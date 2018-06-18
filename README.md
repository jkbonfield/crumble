# Crumble
Exploration of controlled loss of quality values for compressing CRAM files

This directory contains some experimental tools that read a SAM/BAM/CRAM file,
compute which confidence values to keep and which to omit, and emit a new file
with most qualities removed.

See INSTALL for compilation and installation instructions.

**Crumble** uses a simple heterozygous consensus algorithm (taken from
_gap5_) in a couple of modes to produce a consensus call with a
confidence.  If the calls are highly confident then it throws away the
quality values (setting them to fixed high or low depending on whether
they agree with the call), otherwise it keeps them.

To control the output format use -O FMT.  This may also be used to add
additional parameters.  Eg "-O bam,nthreads=8" may be used to enable
multi-threaded BGZF compression when writing the BAM file.

A paper preprint is available with comprehensive benchmarks.

https://www.biorxiv.org/content/early/2018/01/04/243030

