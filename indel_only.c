/*
 * Sets all qualities to the same value Q except for within D bases of
 * any sequence indel.
 */

#define D 20
#define Q 37

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

#include <htslib/sam.h>

//-----------------------------------------------------------------------------
// Main pileup iterator

typedef struct {
    // sam, bam AND hts! All the weasles in a single bag :)
    samFile *fp;
    bam_hdr_t *header;
    hts_itr_t *iter;
} pileup_cd;

int pileup_callback(void *vp, bam1_t *b) {
    pileup_cd *cd = (pileup_cd *)vp;
    return (cd->iter)
	? sam_itr_next(cd->fp, cd->iter, b)
	: sam_read1(cd->fp, cd->header, b);
}

int main(int argc, char **argv) {
    samFile *in, *out = NULL;
    bam_plp_t p_iter;
    int tid, pos;
    int n_plp;
    const bam_pileup1_t *plp;
    bam_hdr_t *header;
    pileup_cd cd;
    hts_itr_t *h_iter = NULL;

    if (argc < 2) {
	fprintf(stderr, "Usage: indel_only SAM/BAM/CRAM-file [region]\n");
	return 1;
    }

    if (!(in = sam_open(argv[1], "rb"))) {
	perror(argv[1]);
	return 1;
    }

    if (!(out = sam_open("-", "w"))) {
	perror("(stdout)");
	return 1;
    }

    if (!(header = sam_hdr_read(in))) {
	fprintf(stderr, "Failed to read file header\n");
	return 1;
    }
    if (out && sam_hdr_write(out, header) != 0) {
	fprintf(stderr, "Failed to write file header\n");
	return 1;
    }

    if (argc > 2) {
        hts_idx_t *idx = sam_index_load(in, argv[1]);
	h_iter = sam_itr_querys(idx, header, argv[2]);
	if (!h_iter || !idx) {
	    fprintf(stderr, "Failed to load index and/or parse iterator.\n");
	    return 1;
	}
	hts_idx_destroy(idx);
    }

    cd.fp = in;
    cd.header = header;
    cd.iter = h_iter;

    p_iter = bam_plp_init(pileup_callback, &cd);
    while ((plp = bam_plp_auto(p_iter, &tid, &pos, &n_plp))) {
	int i;

	if (h_iter) {
	    if (pos < h_iter->beg)
		continue;
	    if (pos >= h_iter->end)
		break;
	}

	for (i = 0; i < n_plp; i++) {
	    bam1_t *b = plp[i].b;

	    if (plp[i].indel || plp[i].is_del) {
		// Mark surrounding D bases
		int x, x_s = plp[i].qpos+1 - D, x_e = plp[i].qpos+1 + D;
		if (x_s < 0) x_s = 0;
		if (x_e >= b->core.l_qseq) x_e = b->core.l_qseq-1;
		for (x = x_s; x <= x_e; x++)
		    bam_get_qual(b)[x] |= 0x80; // marker to keep
	    }
	}

	for (i = 0; i < n_plp; i++) {
	    bam1_t *b = plp[i].b;
	    uint8_t *qual = bam_get_qual(b);
	    if (!plp[i].is_tail) 
		continue;

	    // Correct qualities
	    int x;
	    for (x = 0; x < b->core.l_qseq; x++) {
		if (qual[x] & 0x80)
		    qual[x] &= ~0x80;
		else
		    qual[x] = Q;
	    }
	    sam_write1(out, header, plp[i].b);
	}
    }

    bam_plp_destroy(p_iter);
    bam_hdr_destroy(header);
    if (h_iter) hts_itr_destroy(h_iter);

    if (sam_close(in) != 0) {
	fprintf(stderr, "Error while closing input fd\n");
	return 1;
    }

    if (out && sam_close(out) != 0) {
	fprintf(stderr, "Error while closing output fd\n");
	return 1;
    }

    return 0;
}
