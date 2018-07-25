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

A paper is available with comprehensive benchmarks:

http://dx.doi.org/10.1093/bioinformatics/bty608

Independent benchmarks are also available in a DNAnexus blog post:

https://blog.dnanexus.com/2018-07-23-breaking-down-crumble/
