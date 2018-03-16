#!/bin/sh
#
# Usage compare_vcf A.vcf B.vcf [region-skip.bed]

bcftools=${BCFTOOLS:-bcftools}
vt=${VT:-/nfs/users/nfs_j/jkb/lustre/vt/vt}

if [ $# -lt 2 ]
then
    echo Usage: compare_vcf A.vcf B.vcf [region-skip.bed] [region-include.bed] [region]
    exit 1
fi

v1=$1
v2=$2
exclude=$3
include=$4
region=$5

qual=${QUAL:-30}

href=${HREF:-/nfs/srpipe_references/references/Human/1000Genomes_hs37d5/all/fasta/hs37d5.fa} 

norm() {
    # Also consider norm "-m +both" without a "-d both" step.
    # This produces more differences, but is also more correct.

    v=$1
    if [ x"$region" != x ]
    then
	#$bcftools norm -t $region -f $href $v 2>/dev/null | $bcftools norm -d both -N | $bcftools view -T ^$exclude | $bcftools view -T $include > $v.norm.vcf
	#$bcftools norm -m -both -t $region -f $href $v 2>/dev/null | $bcftools view -T ^$exclude | $bcftools view -T $include > $v.norm.vcf
	if [ x"$exclude" != x ]
	then
	    $bcftools norm -m -both -t $region -f $href $v 2>/dev/null | $vt decompose_blocksub - 2>/dev/null | $bcftools view -T ^$exclude | $bcftools view -T $include > $v.norm.vcf
	else
	    $bcftools norm -m -both -t $region -f $href $v 2>/dev/null | $vt decompose_blocksub - 2>/dev/null | $bcftools view -T $include > $v.norm.vcf
	fi
    elif [ x"$include" != x ]
    then
	$bcftools norm -m -both -f $href $v 2>/dev/null | $vt decompose_blocksub - 2>/dev/null | $bcftools view -T ^$exclude | $bcftools view -T $include > $v.norm.vcf
    elif [ x"$exclude" != x ]
    then
	$bcftools norm -m -both -f $href $v 2>/dev/null | $vt decompose_blocksub - 2>/dev/null | $bcftools view -T ^$exclude > $v.norm.vcf
    else
	$bcftools norm -m -both -f $href $v 2>/dev/null | $vt decompose_blocksub - 2>/dev/null > $v.norm.vcf
    fi
    rm $v.norm.vcf.gz 2>/dev/null
    bgzip $v.norm.vcf
    $bcftools index $v.norm.vcf.gz
}

norm $v1
norm $v2

# NB: consider CHM13_1 chr20:2241427
# The Illumina data is all mapped with high mqual and shows homozygous deletion
# of AAAC in an AAAC STR.  The truth set claims AAAC del is heterozygous.
#
# $Bcftools isec therefore claims it is shared between both sets (in 0002 and 0003.vcf)
# but also only in the query set (0001) due to the extra allele.

#$bcftools isec -p $v1.isec $v1.norm.vcf.gz $v2.norm.vcf.gz
$bcftools isec -c both -p $v2.isec $v1.norm.vcf.gz $v2.norm.vcf.gz
# Or "-c any"? Works better possibly


# 0000 is private to v1 => FN
# 0001 is private to v2 => FP
# 0002/3 are records common to v1/v2 (from v1 or v2 respectively => TP)

# Depth filtering is beneficial to all tools, so we use it in this evaluation.
case $v2 in
    *15x*) DP=30;;
    *)     DP=90;;
esac

case $v2 in
    *gatk*)
	prog=GATK
        #gatk: https://software.broadinstitute.org/gatk/documentation/article.php?id=3225
	# s_filt_exp="QUAL < $qual || QD < 2 || FS > 60 || MQ < 40 || SOR > 3 || MQRankSum < -12.5 || ReadPosRankSum < -8"
	# i_filt_exp="QUAL < $qual || QD < 2 || FS > 200 || ReadPosRankSum < -20"
	# Consider adding DP > 30 here too as it helps a lot.
	s_filt_exp="QUAL < $qual || QD < 2 || FS > 60 || MQ < 40 || SOR > 3 || MQRankSum < -12.5 || ReadPosRankSum < -8 || DP>$DP"
	i_filt_exp="QUAL < $qual || QD < 2 || FS > 200 || ReadPosRankSum < -20 || DP>$DP"
	;;
    *freebayes*)
        prog=Freebayes
        #freebayes: https://wiki.uiowa.edu/download/attachments/145192256/erik%20garrison%20-%20iowa%20talk%202.pdf?api=v2
	s_filt_exp="QUAL < $qual || SAF <= 0 || SAR <= 0 || RPR <= 1 || RPL <= 1 || DP > $DP"
	i_filt_exp=$s_filt_exp
	;;
    *bcftools*)
        prog=Bcftools
        #bcftools: Note QUAL doesn't correlate well for indels.  We also don't have a ratio of depth,
	#so have to manually change this on 15x vs 50x samples here.
	s_filt_exp="QUAL < $qual || DP>$DP"
	i_filt_exp="IDV < 3 || IMF < 0.03 || DP>$DP"
	;;
    *)
        prog=Unknown
	echo Unrecognised type, no specific filters
	s_filt_exp="QUAL<0"
	i_filt_exp="QUAL<0"
	;;
esac

# QUAL 1 is recommended minimum for freebayes to remove detritus.

v1_snp=`    $bcftools view -H -i "TYPE='snp'" $v2.isec/0000.vcf|wc -l`
v2_snp=`    $bcftools view -H -i "TYPE='snp' && QUAL >= 1" $v2.isec/0001.vcf|wc -l`
v2_snp_hq=` $bcftools view -H -i "TYPE='snp' && QUAL >= $qual" $v2.isec/0001.vcf|wc -l`
v2_snp_fi=` $bcftools view    -i "TYPE='snp'" $v2.isec/0001.vcf | bcftools view -H -e "$s_filt_exp" -|wc -l`
v12_snp=`   $bcftools view -H -i "TYPE='snp'" $v2.isec/0002.vcf|wc -l`
v12_snp_hq=`$bcftools view -H -i "TYPE='snp' && QUAL >= $qual" $v2.isec/0003.vcf|wc -l`
v12_snp_fi=`$bcftools view    -i "TYPE='snp'" $v2.isec/0003.vcf | bcftools view -H -e "$s_filt_exp" -|wc -l`

v1_indel=`    $bcftools view -H -i "TYPE='indel'" $v2.isec/0000.vcf|wc -l`
v2_indel=`    $bcftools view -H -i "TYPE='indel' && QUAL >= 1" $v2.isec/0001.vcf|wc -l`
v2_indel_hq=` $bcftools view -H -i "TYPE='indel' && QUAL >= $qual" $v2.isec/0001.vcf|wc -l`
v2_indel_fi=` $bcftools view    -i "TYPE='indel'" $v2.isec/0001.vcf | bcftools view -H -e "$i_filt_exp" -|wc -l`
v12_indel=`   $bcftools view -H -i "TYPE='indel'" $v2.isec/0002.vcf|wc -l`
v12_indel_hq=`$bcftools view -H -i "TYPE='indel' && QUAL >= $qual" $v2.isec/0003.vcf|wc -l`
v12_indel_fi=`$bcftools view    -i "TYPE='indel'" $v2.isec/0003.vcf | bcftools view -H -e "$i_filt_exp" -|wc -l`

# quality trimmed FN aren't the records private to v1 above QUAL, but the
# total number of records not in v12 after filtering.  Thus as we increase
# acceptance threshold to reduce FP we increase FN.
v1_snp_hq=`expr $v1_snp + $v12_snp - $v12_snp_hq`
v1_snp_fi=`expr $v1_snp + $v12_snp - $v12_snp_fi`
v1_indel_hq=`expr $v1_indel + $v12_indel - $v12_indel_hq`
v1_indel_fi=`expr $v1_indel + $v12_indel - $v12_indel_fi`

# Assumption A.vcf is truth set and B.vcf is test set
if [ "$FORMAT" = "tex" ]
then
	printf '\\bigskip\n'
	printf '\\begin{minipage}{\\linewidth}\n'
	printf '\\centering\n'
	printf '\\captionof{table}{%s: FIXME}\n' $prog
	printf '{\\begin{tabular}{ll|r|rr}\n'
	printf 'Variants & & \\textbf{Q>0} & \\textbf{Q>=%d} & \\textbf{Filtered} \\\\ \\midrule\n' $qual
	printf 'SNP   & TP & %7d & %7d & %7d \\\\\n' $v12_snp  $v12_snp_hq  $v12_snp_fi
	printf 'SNP   & FP & %7d & %7d & %7d \\\\\n' $v2_snp   $v2_snp_hq   $v2_snp_fi
	printf 'SNP   & FN & %7d & %7d & %7d \\\\\n' $v1_snp   $v1_snp_hq   $v1_snp_fi
	printf '\\midrule\n';
	printf 'InDel & TP & %7d & %7d & %7d \\\\\n' $v12_indel $v12_indel_hq $v12_indel_fi
	printf 'InDel & FP & %7d & %7d & %7d \\\\\n' $v2_indel  $v2_indel_hq  $v2_indel_fi
	printf 'InDel & FN & %7d & %7d & %7d \\\\\n' $v1_indel  $v1_indel_hq  $v1_indel_fi
	printf '\\end{tabular}}\n'
	printf '\\par\n'
	printf '\\textbf{CRAM qual size x,x}\n'
	printf '\\bigskip\n'
	printf '\\end{minipage}\n'

elif [ "$FORMAT" = "gnuplot" ]
then
	# For gnuplot
	printf "$v2\tall $v12_snp $v2_snp $v1_snp $v12_indel $v2_indel $v1_indel\n"
	printf "$v2\tq$qual $v12_snp_hq $v2_snp_hq $v1_snp_hq $v12_indel_hq $v2_indel_hq $v1_indel_hq\n"
	printf "$v2\tF$qual $v12_snp_fi $v2_snp_fi $v1_snp_fi $v12_indel_fi $v2_indel_fi $v1_indel_fi\n"
else
	printf "SNP          Q>0 /   Q>=$qual / Filtered\n"
	printf "SNP   TP %7d / %7d / %7d\n" $v12_snp  $v12_snp_hq  $v12_snp_fi
	printf "SNP   FP %7d / %7d / %7d\n" $v2_snp   $v2_snp_hq   $v2_snp_fi
	printf "SNP   FN %7d / %7d / %7d\n" $v1_snp   $v1_snp_hq   $v1_snp_fi
	#printf "SNP   %4.1f%% prec, %4.1f%% rec\n" 100.0*$v12_snp_hq/($v12_snp_hq+$v2_snp_hq) 100.0*$v12_snp_hq/($v12_snp_hq+$v1_snp_hq);
	printf "\n";
	printf "InDel TP %7d / %7d / %7d\n" $v12_indel $v12_indel_hq $v12_indel_fi
	printf "InDel FP %7d / %7d / %7d\n" $v2_indel  $v2_indel_hq  $v2_indel_fi
	printf "InDel FN %7d / %7d / %7d\n" $v1_indel  $v1_indel_hq  $v1_indel_fi
	#printf "InDel %4.1f%% prec, %4.1f%% rec\n" 100.0*$v12_indel_hq/($v12_indel_hq+$v2_indel_hq) 100.0*$v12_indel_hq/($v12_indel_hq+$v1_indel_hq);
fi

#rm $v1.norm* $v2.norm*
#rm -rf $v2.isec

