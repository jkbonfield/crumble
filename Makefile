PROGS=crumble

all: $(PROGS)

# From git submodule add -b cram_lossy_names git@github.com:samtools/htslib.git
HTSDIR=./htslib
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a

INCLUDES=-I$(HTSDIR)/include -I$(HTSDIR)
LIBS=-L$(HTSDIR)/lib -L$(HTSDIR) -lhts -Wl,--rpath,$(HTSDIR)/lib -Wl,--rpath,$(HTSDIR) -lpthread -lz -lm -ldl 
#CFLAGS=-g -Wall -Werror
CFLAGS=-O3 -g -Wall -ffast-math -march=native

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

clean:
	-rm *.o $(PROGS)

# Not built by default
indel_only: $(HTSLIB) indel_only.o
	$(CC) -o $@ $< $(CFLAGS) $(LDFLAGS) $(LIBS)

crumble: $(HTSLIB) snp_score.o str_finder.o
	$(CC) -o $@ snp_score.o str_finder.o $(CFLAGS) $(LDFLAGS) $(LIBS)
