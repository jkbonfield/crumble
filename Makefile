PROGS=indel_only snp_score

all: $(PROGS)

#HTSLIB=/software/solexa/pkg/htslib/current/install
HTSLIB=$(HOME)/work/samtools_master/htslib

INCLUDES=-I$(HTSLIB)/include -I$(HTSLIB)
LIBS=-L$(HTSLIB)/lib -L$(HTSLIB) -lhts -Wl,--rpath,$(HTSLIB)/lib -Wl,--rpath,$(HTSLIB) -lm
CFLAGS=-g -Wall

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

clean:
	-rm *.o $(PROGS)

indel_only: indel_only.o
	$(CC) -o $@ $< $(CFLAGS) $(LDFLAGS) $(LIBS)

snp_score: snp_score.o str_finder.o
	$(CC) -o $@ snp_score.o str_finder.o $(CFLAGS) $(LDFLAGS) $(LIBS)
