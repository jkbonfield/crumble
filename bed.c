#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#include <htslib/sam.h>
#include "bed.h"

static int bed_sort(const void *v1, const void *v2) {
    const bed_reg *b1 = (const bed_reg *)v1;
    const bed_reg *b2 = (const bed_reg *)v2;

    if (b1->tid != b2->tid)
	return b1->tid - b2->tid;
    else
	return b1->start - b2->start;
}

// Sort bed regions and Remove duplications / merge overlapping regions.
static void bed_collapse(bed_reg *reg, int *nreg) {
    int nused = *nreg, last_tid = -1, last_end = -1;
    int i, j;

    qsort(reg, nused, sizeof(*reg), bed_sort);

    for (i = j = 0; i < nused; i++) {
	if (reg[i].tid > last_tid ||
	    reg[i].start > last_end) {
	    reg[j++] = reg[i];
	} else if (reg[i].end > reg[j-1].end) {
	    reg[j-1].end = reg[i].end;
	}

	last_tid = reg[i].tid;
	last_end = reg[i].end;
    }
    reg[j++] = reg[i];

    *nreg = j;
}

bed_reg *bed_load(char *fn, bam_hdr_t *header, int *nreg) {
    FILE *fp = fopen(fn, "r");
    bed_reg *reg = NULL;
    int nsz = 0, nused = 0;
    char line[8192], chr[8192];
    int start, end, tid;

    if (!fp) {
	perror(fn);
	goto err;
    }


    while (fgets(line, 8192, fp)) {
	if (strncmp(line, "#", 1) == 0 ||
	    strncmp(line, "track", 5) == 0 ||
	    strncmp(line, "browser", 7) == 0 ||
	    *line == '\n')
	    continue;

	if (sscanf(line, "%s %d %d", chr, &start, &end) != 3) {
	    fprintf(stderr, "Malformed bed line: %s", line);
	    goto err;
	}

	if ((tid = bam_name2id(header, chr)) < 0) {
	    fprintf(stderr, "Unknown reference name: %s\n", chr);
	    goto err;
	}

	if (nused >= nsz) {
	    nsz = nsz ? nsz*2 : 1024;
	    bed_reg *tmp = realloc(reg, nsz * sizeof(*reg));
	    if (!tmp)
		goto err;
	    reg = tmp;
	}

	reg[nused].tid = tid;
	reg[nused].start = start;
	reg[nused].end = end;
	nused++;
    }

    if (fclose(fp) < 0) {
	fp = NULL;
	goto err;
    }

    *nreg = nused;

    //qsort(reg, nused, sizeof(*reg), bed_sort);
    bed_collapse(reg, nreg);

    return reg;

 err:
    free(reg);
    if (fp)
	fclose(fp);
    return NULL;
}

void bed_free(bed_reg *bed) {
    free(bed);
}

#ifdef TEST_MAIN
int main(int argc, char **argv) {
    bed_reg *reg;
    int nreg;

    samFile *in;
    bam_hdr_t *header;

    if (!(in = sam_open(argv[1], "r")))
	return 1;

    if (!(header = sam_hdr_read(in)))
	return 1;

    if (!(reg = bed_load(argv[2], header, &nreg)))
	return 1;

    int i;
    for (i = 0; i < nreg; i++) {
	printf("%d\t%d\t%d\n", reg[i].tid, reg[i].start, reg[i].end);
    }

    sam_close(in);

    return 0;
}
#endif
