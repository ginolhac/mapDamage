import logging
import random

from pathlib import Path

import pysam


_BAM_UNMAPPED = 0x4
_BAM_SECONDARY = 0x100
_BAM_FAILED_QC = 0x200
_BAM_PCR_DUPE = 0x400
_BAM_CHIMERIC = 0x800


class BAMReader:
    def __init__(self, filepath, downsample_to=None, downsample_seed=None):
        self.filepath = Path(filepath)
        self.downsample_to = downsample_to
        self.downsample_seed = downsample_seed

        self.is_stream = False
        if filepath == "-" or self.filepath.is_fifo() or self.filepath.is_char_device():
            self.is_stream = True

        self.handle = pysam.AlignmentFile(self.filepath)

    def close(self):
        self.handle.close()

    def get_references(self):
        return dict(zip(self.handle.references, self.handle.lengths))

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
