/* The MIT License

   Copyright (c) 20082-2012 by Heng Li <lh3@me.com>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <assert.h>
#include <math.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct {
	int n, m;
	uint64_t *a;
} reglist_t;

#include "khash.h"
KHASH_MAP_INIT_STR(reg, reglist_t)
KHASH_SET_INIT_INT64(64)

typedef kh_reg_t reghash_t;

reghash_t *stk_reg_read(const char *fn)
{
	reghash_t *h = kh_init(reg);
	gzFile fp;
	kstream_t *ks;
	int dret;
	kstring_t *str;
	// read the list
	str = calloc(1, sizeof(kstring_t));
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);
	while (ks_getuntil(ks, 0, str, &dret) >= 0) {
		int beg = -1, end = -1;
		reglist_t *p;
		khint_t k = kh_get(reg, h, str->s);
		if (k == kh_end(h)) {
			int ret;
			char *s = strdup(str->s);
			k = kh_put(reg, h, s, &ret);
			memset(&kh_val(h, k), 0, sizeof(reglist_t));
		}
		p = &kh_val(h, k);
		if (dret != '\n') {
			if (ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
				beg = atoi(str->s);
				if (dret != '\n') {
					if (ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
						end = atoi(str->s);
						if (end < 0) end = -1;
					}
				}
			}
		}
		// skip the rest of the line
		if (dret != '\n') while ((dret = ks_getc(ks)) > 0 && dret != '\n');
		if (end < 0 && beg > 0) end = beg, beg = beg - 1; // if there is only one column
		if (beg < 0) beg = 0, end = INT_MAX;
		if (p->n == p->m) {
			p->m = p->m? p->m<<1 : 4;
			p->a = realloc(p->a, p->m * 8);
		}
		p->a[p->n++] = (uint64_t)beg<<32 | end;
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	return h;
}

void stk_reg_destroy(reghash_t *h)
{
	khint_t k;
	if (h == 0) return;
	for (k = 0; k < kh_end(h); ++k) {
		if (kh_exist(h, k)) {
			free(kh_val(h, k).a);
			free((char*)kh_key(h, k));
		}
	}
	kh_destroy(reg, h);
}

/* constant table */

unsigned char seq_nt16_table[256] = {
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15 /*'-'*/,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9,  0,10,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9,  0,10,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

unsigned char seq_nt6_table[256] = {
    0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,

    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

char *seq_nt16_rev_table = "XACMGRSVTWYHKDBN";
unsigned char seq_nt16to4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };
unsigned char seq_nt16comp_table[] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };
int bitcnt_table[] = { 4, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };
char comp_tab[] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};



/* composition */
int stk_comp(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int l, c, upper_only = 0;
	reghash_t *h = 0;
	reglist_t dummy;

	while ((c = getopt(argc, argv, "ur:")) >= 0) {
		switch (c) {
			case 'u': upper_only = 1; break;
			case 'r': h = stk_reg_read(optarg); break;
		}
	}
	if (argc == optind && isatty(fileno(stdin))) {
		fprintf(stderr, "Usage:  seqtk comp [-u] [-r in.bed] <in.fa>\n\n");
		fprintf(stderr, "Output format: chr, length, #A, #C, #G, #T, #2, #3, #4, #CpG, #tv, #ts, #CpG-ts\n");
		return 1;
	}
	fp = optind < argc && strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	dummy.n= dummy.m = 1; dummy.a = calloc(1, 8);
	while ((l = kseq_read(seq)) >= 0) {
		int i, k;
		reglist_t *p = 0;
		if (h) {
			khint_t k = kh_get(reg, h, seq->name.s);
			if (k != kh_end(h)) p = &kh_val(h, k);
		} else {
			p = &dummy;
			dummy.a[0] = l;
		}
		for (k = 0; p && k < p->n; ++k) {
			int beg = p->a[k]>>32, end = p->a[k]&0xffffffff;
			int la, lb, lc, na, nb, nc, cnt[11];
			if (beg > 0) la = seq->seq.s[beg-1], lb = seq_nt16_table[la], lc = bitcnt_table[lb];
			else la = 'a', lb = -1, lc = 0;
			na = seq->seq.s[beg]; nb = seq_nt16_table[na]; nc = bitcnt_table[nb];
			memset(cnt, 0, 11 * sizeof(int));
			for (i = beg; i < end; ++i) {
				int is_CpG = 0, a, b, c;
				a = na; b = nb; c = nc;
				na = seq->seq.s[i+1]; nb = seq_nt16_table[na]; nc = bitcnt_table[nb];
				if (b == 2 || b == 10) { // C or Y
					if (nb == 4 || nb == 5) is_CpG = 1;
				} else if (b == 4 || b == 5) { // G or R
					if (lb == 2 || lb == 10) is_CpG = 1;
				}
				if (upper_only == 0 || isupper(a)) {
					if (c > 1) ++cnt[c+2];
					if (c == 1) ++cnt[seq_nt16to4_table[b]];
					if (b == 10 || b == 5) ++cnt[9];
					else if (c == 2) {
						++cnt[8];
					}
					if (is_CpG) {
						++cnt[7];
						if (b == 10 || b == 5) ++cnt[10];
					}
				}
				la = a; lb = b; lc = c;
			}
			if (h) printf("%s\t%d\t%d", seq->name.s, beg, end);
			else printf("%s\t%d", seq->name.s, l);
			for (i = 0; i < 11; ++i) printf("\t%d", cnt[i]);
			putchar('\n');
		}
		fflush(stdout);
	}
	free(dummy.a);
	kseq_destroy(seq);
	stk_reg_destroy(h);
	gzclose(fp);
	return 0;
}

static void print_seq(FILE *fpout, const kseq_t *ks, int begin, int end)
{
	int i;
	if (begin >= end) return; // FIXME: why may this happen? Understand it!
	fprintf(fpout, "%c%s:%d-%d", ks->qual.l? '@' : '>', ks->name.s, begin+1, end);
	for (i = begin; i < end && i < ks->seq.l; ++i) {
		if ((i - begin)%60 == 0) fputc('\n', fpout);
		fputc(ks->seq.s[i], fpout);
	}
	fputc('\n', fpout);
	if (ks->qual.l == 0) return;
	fputs("+\n", fpout);
	for (i = begin; i < end && i < ks->qual.l; ++i) {
		if ((i - begin)%60 == 0) fputc('\n', fpout);
		fputc(ks->qual.s[i], fpout);
	}
	fputc('\n', fpout);
}

/* main function */
static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   seqtk <command> <arguments>\n");
	fprintf(stderr, "Version: 1.0-r82-dirty\n\n");
	fprintf(stderr, "Command: comp      get the nucleotide composition of FASTA/Q\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc == 1) return usage();
	if (strcmp(argv[1], "comp") == 0) stk_comp(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'. Abort!\n", argv[1]);
		return 1;
	}
	return 0;
}
