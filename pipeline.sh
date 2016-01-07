#!/bin/sh

# Tests the effect of lossy quality compression.
# Usage:  pipeline foo.bam ref.fa truth.vcf.gz

bcftools=bcftools
samtools=samtools
bgzip=bgzip
freebayes=/software/solexa/pkg/freebayes/1.0.1/bin/freebayes

if [ $# -ne 3 ]
then
    echo "Usage: pipeline foo.bam ref.fa truth.vcf.gz"
    exit 1
fi

bam=$1
ref=$2
truth=$3

prefix=__$$
echo prefix: $prefix

# $bcftools index $truth

#-----------------------------------------------------------------------------
# With freebayes; better on indels
$samtools view -u $bam | $freebayes -f $ref - | bcftools norm -f $ref | $bgzip > $prefix.freebayes.norm.vcf.gz
$bcftools index $prefix.freebayes.norm.vcf.gz
$bcftools stats -s- $truth $prefix.freebayes.norm.vcf.gz > $prefix.freebayes.stats

#-----------------------------------------------------------------------------
# With samtools; better on SNPs

#$samtools mpileup -g -f $ref $bam | $bcftools call -vmO z -o $prefix.vcf.gz
$samtools mpileup -m 2 -p -F 0.1 -g -f $ref $bam | $bcftools call -vmO z -o $prefix.vcf.gz
$bcftools norm -f $ref $prefix.vcf.gz | $bgzip > $prefix.samtools.norm.vcf.gz
$bcftools index $prefix.samtools.norm.vcf.gz
$bcftools stats -s- $truth $prefix.samtools.norm.vcf.gz > $prefix.samtools.stats

#-----------------------------------------------------------------------------
# With GATK
