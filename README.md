# Crumble
Exploration of controlled loss of quality values for compressing CRAM files

This directory contains some experimental tools that read a SAM/BAM/CRAM file,
compute which confidence values to keep and which to omit, and emit a new file
with most qualities removed.

It depends on a pre-compiled htslib.  Preferably this will be against
an installed copy, in which case you need to specify where to find it
via e.g.:

    make CPPFLAGS="-I/usr/local/htslib/include" LDFLAGS="-L/usr/local/htslib/lib"

As a simpler alternative, it will default to using a sibling ../htslib
source tree, but this still needs compiling first.

**Crumble** uses a simple heterozygous consensus algorithm (taken from
_gap5_) in a couple of modes to produce a consensus call with a
confidence.  If the calls are highly confident then it throws away the
quality values (setting them to fixed high or low depending on whether
they agree with the call), otherwise it keeps them.

It also keeps all quality values for reads with mapping quality 0 (they may need
to be realigned elsewhere) and within the proximity of any uncertain indel calls,
where "proximity" is determined by a simply short tandem repeat search.

To control the output format use -O FMT.  This may also be used to add
additional parameters.  Eg "-O bam,nthreads=8" may be used to enable
multi-threaded BGZF compression when writing the BAM file.

**pipeline.sh** is a noddy shell script to run samtools mpileup / bcftools call
and freebayes on a sam/bam file and compares the results against a known truth set.
Earlier results from some of these can be seen in BENCHMARKS_chr20, although this is
now a little out of date.


