#!/bin/sh

# Tests the effect of lossy quality compression.
# Usage:  pipeline foo.bam ref.fa truth.vcf.gz

bcftools=bcftools
#bcftools=bcftools_pd3
samtools=samtools
bgzip=bgzip
freebayes=/software/solexa/pkg/freebayes/1.0.1/bin/freebayes
gatk_jar=/software/solexa/pkg/GATK/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar
gatk_key=/software/vertres/installs/gatk/vertres_gatk.key

if [ $# -ne 3 ]
then
    echo "Usage: pipeline foo.bam ref.fa truth.vcf.gz"
    exit 1
fi

bam=$1
ref=$2
truth=$3

prefix=$bam
echo prefix: $prefix

# $bcftools index $truth

#-----------------------------------------------------------------------------
# With freebayes; better on indels
($samtools view -u $bam | $freebayes -jH0 -f $ref - | bcftools norm -f $ref | $bgzip > $prefix.freebayes.norm.vcf.gz
$bcftools index $prefix.freebayes.norm.vcf.gz
$bcftools stats -s- $truth $prefix.freebayes.norm.vcf.gz > $prefix.freebayes.stats) &

#-----------------------------------------------------------------------------
# With samtools; better on SNPs

($samtools mpileup -m 2 -p -F 0.1 -gu -f $ref $bam | $bcftools call -vmO v | $bcftools norm -f $ref | bgzip > $prefix.samtools.norm.vcf.gz
#$bcftools mpileup -m 2 -p -F 0.1 -gu -f $ref $bam | $bcftools call -vmO z -o $prefix.vcf.gz
#$bcftools norm -f $ref $prefix.vcf.gz | $bgzip > $prefix.samtools.norm.vcf.gz
$bcftools index $prefix.samtools.norm.vcf.gz
$bcftools stats -s- $truth $prefix.samtools.norm.vcf.gz > $prefix.samtools.stats) &

#-----------------------------------------------------------------------------
# With GATK
samtools index $bam
(/software/bin/java -Xmx4g -jar $gatk_jar -T HaplotypeCaller -R $ref -I $bam --genotyping_mode DISCOVERY -stand_emit_conf 0 -stand_call_conf 30 -K $gatk_key | bcftools norm -f $ref - | bgzip > $prefix.gatk.norm.vcf.gz
$bcftools index $prefix.gatk.norm.vcf.gz
$bcftools stats -s- $truth $prefix.gatk.norm.vcf.gz > $prefix.gatk.stats) &

wait

