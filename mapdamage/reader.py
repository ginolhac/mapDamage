import logging
import random

from pathlib import Path

import pysam


_BAM_UNMAPPED = 0x4
_BAM_SECONDARY = 0x100
_BAM_FAILED_QC = 0x200
_BAM_PCR_DUPE = 0x400
_BAM_CHIMERIC = 0x800


class BAMError(RuntimeError):
    pass


class BAMReader:
    def __init__(
        self, filepath, merge_libraries=False, downsample_to=None, downsample_seed=None,
    ):
        log = logging.getLogger(__name__)

        self.filepath = Path(filepath)
        self.downsample_to = downsample_to
        self.downsample_seed = downsample_seed

        self.is_stream = False
        if filepath == "-" or self.filepath.is_fifo() or self.filepath.is_char_device():
            self.is_stream = True

        self.handle = pysam.AlignmentFile(self.filepath)

        self._merge_libraries = merge_libraries
        self._readgroups = {}
        self._libraries = {}

        if merge_libraries:
            self._readgroups[None] = ("*", "*")
            self._libraries[("*", "*")] = set((None,))
        else:
            self._readgroups = self._collect_readgroups(log, self.handle)
            for readgroup, library in self._readgroups.items():
                self._libraries.setdefault(library, set()).add(readgroup)

        log.info("Found %i libraries in BAM file", len(self._libraries))

    def close(self):
        self.handle.close()

    def get_references(self):
        return dict(zip(self.handle.references, self.handle.lengths))

    def get_libraries(self):
        return self._libraries.keys()

    def get_sample_and_library(self, read):
        if self._merge_libraries:
            return self._readgroups[None]

        try:
            readgroup = read.get_tag("RG")
        except KeyError:
            raise BAMError(
                "Read %r has no read-group. Either fix BAM or use --merge-libraries"
                % (read.query_name,)
            )

        try:
            return self._readgroups[readgroup]
        except KeyError:
            raise BAMError(
                "Read %r has read-group not listed in BAM header (%r); either fix BAM "
                "or use --merge-libraries" % (read.query_name, readgroup)
            )

    def __iter__(self):
        log = logging.getLogger(__name__)
        if self.downsample_to is None:
            return self._filter_reads(self.handle)
        elif self.downsample_to < 1:
            log.debug("Downsampling BAM to %.1f%%", self.downsample_to * 100)
            return self._downsample_to_fraction(
                self.handle, self.downsample_to, self.downsample_seed
            )
        else:
            log.debug("Downsampling BAM to %d random reads", self.downsample_to)
            return self._downsample_to_fixed_number(
                self.handle, self.downsample_to, self.downsample_seed
            )

    @classmethod
    def _collect_readgroups(cls, log, handle):
        readgroups = {}
        for readgroup in handle.header.get("RG", ()):
            log.debug(
                "Found readgroup %r with SM=%r and LB=%r",
                readgroup.get("ID"),
                readgroup.get("SM"),
                readgroup.get("LB"),
            )

            try:
                readgroups[readgroup["ID"]] = (readgroup["SM"], readgroup["LB"])
            except KeyError as error:
                raise BAMError(
                    "Incomplete readgroup found: %s is missing %s. "
                    "Either fix BAM or use --merge-libraries"
                    % (readgroup.get("ID", "Unnamed readgroup"), error)
                )

        return readgroups

    @classmethod
    def _filter_reads(cls, handle):
        filtered_flags = (
            _BAM_UNMAPPED
            | _BAM_SECONDARY
            | _BAM_FAILED_QC
            | _BAM_PCR_DUPE
            | _BAM_CHIMERIC
        )

        for read in handle:
            if not (read.flag & filtered_flags):
                yield read

    @classmethod
    def _downsample_to_fraction(cls, handle, downsample_to, seed):
        if not (0 <= downsample_to < 1):
            raise ValueError(downsample_to)

        rand = random.Random(seed)
        for read in cls._filter_reads(handle):
            if rand.random() < downsample_to:
                yield read

    @classmethod
    def _downsample_to_fixed_number(cls, handle, downsample_to, seed):
        # use reservoir sampling
        if downsample_to < 1:
            raise ValueError(downsample_to)

        downsample_to = int(downsample_to)
        rand = random.Random(seed)
        sample = [None] * downsample_to
        for (index, record) in enumerate(cls._filter_reads(handle)):
            if index >= downsample_to:
                index = rand.randint(0, index)
                if index >= downsample_to:
                    continue
            sample[index] = record

        result = [read for read in sample if read is not None]
        # Sampling does not preserve input order
        result.sort(key=lambda read: (read.reference_id, read.reference_start))

        return iter(result)
