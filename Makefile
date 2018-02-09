PROGS=crumble

all: $(PROGS)

LIBS ?= -lhts -lpthread -lz -lm -ldl
#CFLAGS?=-g -Wall -Werror
CFLAGS?=-O3 -g -Wall -ffast-math -march=native

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $<

clean:
	-rm *.o $(PROGS)

# Not built by default
indel_only: indel_only.o
	$(CC) -o $@ $< $(CFLAGS) $(LDFLAGS) $(LIBS)

crumble: snp_score.o str_finder.o
	$(CC) -o $@ snp_score.o str_finder.o $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LIBS)
