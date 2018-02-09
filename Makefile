PROGS=crumble
PREFIX=/usr/local

all: $(PROGS)

#CFLAGS=-g -Wall -Werror
CFLAGS=-O3 -g -Wall -ffast-math -march=native
LIBS=-lhts -lpthread -lz -lm -ldl

# Defaults for building against htslib in a sibling directory.
# These have no include or lib suffixes as we are using the source tree.
# (This needs to have been compiled first.)
CPPFLAGS=-I../htslib
LDFLAGS=-L../htslib -Wl,--rpath,../htslib

# To compile against an installed (non-source) htslib, use e.g.
# make CPPFLAGS="-I/usr/local/htslib/include" LDFLAGS="-L/usr/local/htslib/lib -Wl,--rpath,/usr/local/htslib/lib"

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $<

clean:
	-rm *.o $(PROGS)

# Not built by default
indel_only: indel_only.o
	$(CC) -o $@ $< $(CFLAGS) $(LDFLAGS) $(LIBS)

crumble: snp_score.o str_finder.o
	$(CC) -o $@ snp_score.o str_finder.o $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LIBS)

install: $(PROGS)
	cp $(PROGS) $(PREFIX)/bin
