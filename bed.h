#ifndef BED_H
#define BED_H

typedef struct reg {
    int tid , start, end;
} bed_reg;

bed_reg *bed_load(char *fn, bam_hdr_t *header, int *nreg);
void bed_free(bed_reg *bed);

#endif /* BED_LOAD_H */
