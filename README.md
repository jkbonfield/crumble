# Crumble
Exploration of controlled loss of quality values for compressing CRAM files

This directory contains some experimental tools that read a SAM/BAM/CRAM file,
compute which confidence values to keep and which to omit, and emit a new file
with most qualities removed.

**Indel_only** is a noddy exploration that simply removes all qualities except
those within _D_ distance of an indel.

**Crumble** performs a more complete trial.  It uses a simple heterozygous
consensus algorithm (ripped out of _gap5_) in a couple of modes to produce a
consensus call with a confidence.  If the calls are highly confident then it
throws away the quality values (setting them to fixed high or low depending on
whether they agree with the call), otherwise it keeps them.

It also keeps all quality values for reads with mapping quality 0 (they may need
to be realigned elsewhere) and within the proximity of any uncertain indel calls,
where "proximity" is determined by a simply short tandem repeat search.

**pipeline.sh** is a noddy shell script to run samtools mpileup / bcftools call
and freebayes on a sam/bam file and compares the results against a known truth set.
Earlier results from some of these can be seen in BENCHMARKS_chr20, although this is
now a little out of date.


