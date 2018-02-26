// Various extra heuristics for storing qualities.  The general rule
// of thumb is that if we think the alignment is suspect in some way
// then we need to keep all the qualities so that it can be remapped
// to a new place.  Consider migration from GRCh37 to h38 - even if it
// has high mapping quality in GRCh37, if we can detect something
// suspicious then we preserve our accuracy when we subsequently
// remap.
//
// Additionally we may be able to generate a bed file showing regions
// that shouldn't be considered as high quality for snp/indel
// calling.


// TODO: excessive depth => mismapped. (keep quals) Eg
// 16:46393052

// TODO: More than 2 haplotypes for 1 single sample. 
// 1:142537663 (Also extra deep)
// 1:8404952
// 17:22251551

// TODO: Prevelance of lots of low mapping quality may mean other
// high quality reads should also be considered low mapping quality?
// (But perhaps only if they're not softclipped.)  Lack of
// soft-clipping implies collapsed repeat, so casts doubt on all data.
// DONE: need to test
// 1:8404952

// TODO: Concordant soft-clip detection
// Regions with lots of soft-clip that aligns well against each other
// indicate a large insertion or possibly an expanded repeat.  These
// sometimes give rise to fake SNPs.
// [Less complex solution, count %age of reads with clip-point larger
// than a specific length at a specific site. If > X then preserve all
// reads with a clip point at that site (regardless of length).  Or
// all reads, not just those with soft-clip?]
// 1:231337816
// 18:312676
// X:104252683

// TODO: Insertion concordance.
// This is hard to do with samtools pileup as we only get ref by ref
// rather than consensus by consensus, but if there is a lot of
// variance between the inserted bases then it implies we have suspect
// alignments.
// DONE: Trivial implementation is simply length concordance; if there
// are many different lengths then it's probably suspect.

//#define DEBUG

#define CRUMBLE_VERSION "0.8"

/*
 * Prunes quality based on snp calling score.
 * 
 * Bi-allelic consensus is computed (possibly as hom).
 * If call has high confidence, then bin to 2 qualities.
 * - high if base is one of the 1-2 consensus calls.
 * - low otherwise.
 *
 * Exceptions.
 * - Within 'D' bases of any read indel.
 * - Within 'D' bases of a soft-clip (possible clipped indel).
 * - Any sequence with mapping quality of <= M.
 * - Any column with high discrepancy, implying possibly tri-allelic.
 *
 */

/*
 * TODO: Het indels may mean we quantise the indel but not the
 * non-indel.  This gives a bias.  Need to treat both the same.
 *
 * Identify low complexity regions.  If we're preserving confidence
 * then it needs to be for all bases in that low-complexity section so
 * that local realignment doesn't change qualities.
 *
 * Consider using STR method for detecting start/end range of bases to
 * keep for SNPs as well as indels, particularly when the sequence is
 * soft-clipped during an STR.  (Possibility of misaligned bases
 * then.)  Or just a low complexity filter? Or concordant soft-clips?
 *
 * For indels, consider the soft-clip adjustment on reads ending in
 * STR adjacent to indels as these bases don't confirm the count.
 */

// Default params
#define MAX_DEPTH 20000

#define QL 5
#define QM 25 // below => QL, else QH
#define QH 40

// If mapping qual <= MIN_MQUAL we preserve all quality values.
#define MIN_MQUAL 0

// Whether to allow quality reduction on mismatching bases (ie QL)
#define REDUCE_QUAL 1
#define BINARY_QUAL 0

// Standard gap5 algorithm; set MIN_QUAL_A to 0 to disable
#define MIN_QUAL_A 0
#define MIN_INDEL_A 50
#define MIN_DISCREP_A 2.0

// With mqual adjustment; set MIN_QUAL_B to 0 to disable
#define MIN_QUAL_B 75
#define MIN_INDEL_B 150
#define MIN_DISCREP_B 1.0

//#define MIN_QUAL_B 50
//#define MIN_INDEL_B 100
//#define MIN_DISCREP_B 1.5

// Extra growth to expand indel qual region.
// New region = old_region +/- (region_len*STR_MUL + STR_ADD)
#define I_STR_MUL 1.1
#define S_STR_MUL 0.0

#define I_STR_ADD 2
#define S_STR_ADD 0

// Prevalence of low mapping quality, > PERC => store all
// Lower => larger files
#define LOW_MQUAL_PERC 0.5

// Amount of variable sized insertion
#define INS_LEN_PERC 0.1

// Percentage of seqs with large soft-clips
#define CLIP_PERC 0.2

// Amount of over-depth to consider this as as suspect region
#define OVER_DEPTH 3.0

// Percentage of reads spanning indel.
#define INDEL_OVERLAP_PERC 0.5

#define BED_DIST 50

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <getopt.h>
#include <inttypes.h>
#include <unistd.h>

#include <htslib/sam.h>
#include <htslib/cram.h>
#include <htslib/khash.h>

#include "str_finder.h"
#include "tree.h"

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

#ifndef ABS
#define ABS(a) ((a)>=0?(a):-(a))
#endif

KHASH_SET_INIT_INT(aux_exists)
typedef khash_t(aux_exists) *auxhash_t;

typedef struct {
    int    reduce_qual, binary_qual;
    int    iSTR_add,  sSTR_add;
    double iSTR_mul, sSTR_mul;
    int    qlow, qcutoff, qhigh;
    int    min_mqual;
    char  *region;

    // Standard gap5 algorithm
    int    min_qual_A;
    int    min_indel_A;
    double min_discrep_A;

    // Mqual adjusted algorithm
    int    min_qual_B;
    int    min_indel_B;
    double min_discrep_B;

    // Tag white/black lists
    auxhash_t aux_whitelist;
    auxhash_t aux_blacklist;

    double low_mqual_perc;
    double clip_perc;
    double ins_len_perc;
    double over_depth;
    double indel_ov_perc;
    FILE *bed_fp;
    int verbose;
    int pblock;
    int softclip;
    int noPG;

    // For BD/BI tag adjustments
    int BD_low, BD_mid, BD_high;
    int BI_low, BI_mid, BI_high;
} cram_lossy_params;

//-----------------------------------------------------------------------------
// Binning

static int bin2[256];

void init_bins(cram_lossy_params *p) {
    int i;

    // Binary mode - just high/low
    for (i = 0; i < p->qcutoff; i++)
	bin2[i] = p->qlow;
    for (; i < 256; i++)
	bin2[i] = p->qhigh;
}

//-----------------------------------------------------------------------------
// Bits ripped out of gap5's consensus algorithm; an interim solution

#define CONS_DISCREP    4
#define CONS_ALL        15

#define CONS_MQUAL      16

typedef struct {
    /* the most likely base call - we never call N here */
    /* A=0, C=1, G=2, T=3, *=4 */
    int call;

    /* The most likely heterozygous base call */
    /* Use "ACGT*"[het / 5] vs "ACGT*"[het % 5] for the combination */
    int het_call;

    /* Log-odds for het_call */
    int het_phred;

    /* Single phred style call */
    unsigned char phred;

    /* Sequence depth */
    int depth;

    /* Discrepancy search score */
    float discrep;
} consensus_t;

#define P_HET 1e-6

#define LOG10            2.30258509299404568401
#define TENOVERLOG10     4.34294481903251827652
#define TENLOG2OVERLOG10 3.0103


/* Sequencing technologies for seq_t.seq_tech; 5 bits, so max=31 */
#define STECH_UNKNOWN    0
#define STECH_SANGER     1
#define STECH_SOLEXA     2
#define STECH_SOLID      3
#define STECH_454        4
#define STECH_HELICOS    5
#define STECH_IONTORRENT 6
#define STECH_PACBIO     7
#define STECH_ONT        8
#define STECH_LAST       8 // highest value

double tech_undercall[] = {
    1.00, // unknown
    1.00, // sanger
    1.00, // solexa/illumina
    1.00, // solid
    1.00, // 454
    1.00, // helicos
    1.00, // iontorrent
    1.00, // pacbio
    1.63, // ont
};

#ifdef __GNUC__
#define ALIGNED(x) __attribute((aligned(x)))
#else
#define ALIGNED(x)
#endif

static double prior[25]    ALIGNED(16);  /* Sum to 1.0 */
static double lprior15[15] ALIGNED(16);  /* 15 combinations of {ACGT*} */

/* Precomputed matrices for the consensus algorithm */
static double pMM[9][101] ALIGNED(16);
static double p__[9][101] ALIGNED(16);
static double p_M[9][101] ALIGNED(16);
static double po_[9][101] ALIGNED(16);
static double poM[9][101] ALIGNED(16);
static double poo[9][101] ALIGNED(16);
static double puu[9][101] ALIGNED(16);
static double pum[9][101] ALIGNED(16);
static double pmm[9][101] ALIGNED(16);

static double e_tab_a[1002]  ALIGNED(16);
static double *e_tab = &e_tab_a[500];
static double e_tab2_a[1002] ALIGNED(16);
static double *e_tab2 = &e_tab2_a[500];
static double e_log[501]     ALIGNED(16);

/*
 * Lots of confusing matrix terms here, so some definitions will help.
 *
 * M = match base
 * m = match pad
 * _ = mismatch
 * o = overcall
 * u = undercall
 *
 * We need to distinguish between homozygous columns and heterozygous columns,
 * done using a flat prior.  This is implemented by treating every observation
 * as coming from one of two alleles, giving us a 2D matrix of possibilities
 * (the hypotheses) for each and every call (the observation).
 *
 * So pMM[] is the chance that given a call 'x' that it came from the
 * x/x allele combination.  Similarly p_o[] is the chance that call
 * 'x' came from a mismatch (non-x) / overcall (consensus=*) combination.
 *
 * Examples with observation (call) C and * follows
 *
 *  C | A  C  G  T  *          * | A  C  G  T  * 
 *  -----------------	       ----------------- 
 *  A | __ _M __ __ o_	       A | uu uu uu uu um
 *  C | _M MM _M _M oM	       C | uu uu uu uu um
 *  G | __ _M __ __ o_	       G | uu uu uu uu um
 *  T | __ _M __ __ o_	       T | uu uu uu uu um
 *  * | o_ oM o_ o_ oo	       * | um um um um mm
 *
 * In calculation terms, the _M is half __ and half MM, similarly o_ and um.
 *
 * Relative weights of substitution vs overcall vs undercall are governed on a
 * per base basis using the P_OVER and P_UNDER scores (subst is 1-P_OVER-P_UNDER).
 *
 * The heterozygosity weight though is a per column calculation as we're
 * trying to model whether the column is pure or mixed. Hence this is done
 * once via a prior and has no affect on the individual matrix cells.
 */

static void consensus_init(double p_het) {
    int i, t;

    for (i = -500; i <= 500; i++)
    	e_tab[i] = exp(i);
    for (i = -500; i <= 500; i++)
    	e_tab2[i] = exp(i/10.);
    for (i = 0; i <= 500; i++)
	e_log[i] = log(i);

    // Heterozygous locations
    for (i = 0; i < 25; i++)
	prior[i] = p_het / 20;
    prior[0] = prior[6] = prior[12] = prior[18] = prior[24] = (1-p_het)/5;

    lprior15[0]  = log(prior[0]);
    lprior15[1]  = log(prior[1]*2);
    lprior15[2]  = log(prior[2]*2);
    lprior15[3]  = log(prior[3]*2);
    lprior15[4]  = log(prior[4]*2);
    lprior15[5]  = log(prior[6]);
    lprior15[6]  = log(prior[7]*2);
    lprior15[7]  = log(prior[8]*2);
    lprior15[8]  = log(prior[9]*2);
    lprior15[9]  = log(prior[12]);
    lprior15[10] = log(prior[13]*2);
    lprior15[11] = log(prior[14]*2);
    lprior15[12] = log(prior[18]);
    lprior15[13] = log(prior[19]*2);
    lprior15[14] = log(prior[24]);


    // Rewrite as new form
    for (t = STECH_UNKNOWN; t <= STECH_LAST; t++) {
	for (i = 1; i < 101; i++) {
	    double prob = 1 - pow(10, -i / 10.0);

//	    if (t == STECH_ONT)
//		prob = 0.85; // Fake FIXED prob for now

	    // May want to multiply all these by 5 so pMM[i] becomes close
	    // to -0 for most data. This makes the sums increment very slowly,
	    // keeping bit precision in the accumulator.
#if 0
	    double p_overcall[] = {
		0.001, // unknown
		0.010, // sanger
		0.001, // solexa/illumina
		0.001, // solid
		0.010, // 454
		0.010, // helicos
		0.010, // iontorrent
		0.010, // pacbio
		0.050, // ont
	    };

	    double p_undercall[] = {
		0.001, // unknown
		0.010, // sanger
		0.001, // solexa/illumina
		0.001, // solid
		0.010, // 454
		0.010, // helicos
		0.010, // iontorrent
		0.010, // pacbio
		0.280, // ont
	    };

	    double norm = (1-p_overcall[t])*prob + 3*((1-p_overcall[t])*(1-prob)/3)
		+ p_overcall[t]*(1-prob);
	    pMM[t][i] = log((1-p_overcall[t]) * prob /norm);
	    p__[t][i] = log((1-p_overcall[t]) * (1-prob)/3 /norm);
	    poo[t][i] = log((p_overcall[t]*(1-prob)) /norm);

	    p_M[t][i] = log((exp(pMM[t][i]) + exp(p__[t][i]))/2);
	    po_[t][i] = log((exp(p__[t][i]) + exp(poo[t][i]))/2);
	    poM[t][i] = log((exp(pMM[t][i]) + exp(poo[t][i]))/2);

	    // [t]* observation vs base
	    norm = p_undercall[t]*(1-prob)*4 + (1-p_undercall[t])*prob;
	    puu[t][i] = log((p_undercall[t] * (1-prob)) /norm);
	    pmm[t][i] = log((1-p_undercall[t])*prob /norm);
	    pum[t][i] = log((exp(puu[t][i]) + exp(pmm[t][i]))/2);
#else
	    //prob = 1-(1-prob)*2; if (prob < 0.1) prob = 0.1; // Fudge

	    pMM[t][i] = log(prob/5);
	    p__[t][i] = log((1-prob)/20);
	    p_M[t][i] = log((exp(pMM[t][i]) + exp(p__[t][i]))/2);

	    puu[t][i] = p__[t][i];

	    poM[t][i] = p_M[t][i] *= tech_undercall[t];
	    po_[t][i] = p__[t][i] *= tech_undercall[t];
	    poo[t][i] = p__[t][i] *= tech_undercall[t];
	    pum[t][i] = p_M[t][i] *= tech_undercall[t];
	    pmm[t][i] = pMM[t][i] *= tech_undercall[t];
#endif
	}

	pMM[t][0] = pMM[t][1];
	p__[t][0] = p__[t][1];
	p_M[t][0] = p_M[t][1];

	pmm[t][0] = pmm[t][1];
	poo[t][0] = poo[t][1];
	po_[t][0] = po_[t][1];
	poM[t][0] = poM[t][1];
	puu[t][0] = puu[t][1];
	pum[t][0] = pum[t][1];
    }
}

static inline double fast_exp(double y) {
    if (y >= -50 && y <= 50)
	return e_tab2[(int)(y*10)];

    if (y < -500)
	y = -500;
    if (y > 500)
	y = 500;

    //    printf("%f => %g %g\n", y, exp(y), e_tab[(int)y]);

    return e_tab[(int)y];
}

/*Taylor (deg 3) implementation of the log: http://www.flipcode.com/cgi-bin/fcarticles.cgi?show=63828*/
inline double fast_log2 (double val)
{
   register int64_t *const     exp_ptr = ((int64_t*)&val);
   register int64_t            x = *exp_ptr;
   register const int      log_2 = ((x >> 52) & 2047) - 1024;
   x &= ~(2047LL << 52);
   x += 1023LL << 52;
   *exp_ptr = x;

   val = ((-1.0f/3) * val + 2) * val - 2.0f/3;

   return val + log_2;
}

inline double fast_log (double val) {
    return fast_log2(val)*0.69314718;
}

//#define fast_exp exp
//#define fast_log log

#define ph_log(x) (-TENLOG2OVERLOG10*fast_log2((x)))


/*
 * As per calculate_consensus_bit_het but for a single pileup column.
 */
int calculate_consensus_pileup(int flags,
			       const bam_pileup1_t *p,
			       int np,
			       consensus_t *cons) {
    int i, j;
    static int init_done =0;
    static double q2p[101], mqual_pow[256];
    double min_e_exp = DBL_MIN_EXP * log(2) + 1;

    double S[15] ALIGNED(16) = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double sumsC[6] = {0,0,0,0,0,0}, sumsE = 0;
    int depth = 0;

    /* Map the 15 possible combinations to 1-base or 2-base encodings */
    static int map_sing[15] ALIGNED(16) =
	{0, 5, 5, 5, 5,
	    1, 5, 5, 5,
	       2, 5, 5,
	          3, 5,
	             4};
    static int map_het[15] ALIGNED(16) =
	{0,  1,  2,  3,  4,
	     6,  7,  8,  9,
	        12, 13, 14,
	            18, 19,
	                24};

    if (!init_done) {
	init_done = 1;
	consensus_init(P_HET);

	for (i = 0; i <= 100; i++) {
	    q2p[i] = pow(10, -i/10.0);
	}

	for (i = 0; i < 255; i++) {
	    //mqual_pow[i] = 1-pow(10, -(i+.01)/10.0);
	    //mqual_pow[i] = 1-pow(10, -(i/3+.1)/10.0);
	    mqual_pow[i] = 1-pow(10, -(i/2+.05)/10.0);
	}
	// unknown mqual
	mqual_pow[255] = mqual_pow[10];
    }

    /* Initialise */
    int counts[6] = {0};

    /* Accumulate */
    int n;
    //printf("-----\n");

    // FIXME: also seed with unknown alleles so low coverage data is
    // less confident.

    for (n = 0; n < np; n++) {
	bam1_t *b = p[n].b;
	uint8_t base = bam_seqi(bam_get_seq(b), p[n].qpos);
	uint8_t qual = bam_get_qual(b)[p[n].qpos];
	const int stech = STECH_SOLEXA;

	// =ACM GRSV TWYH KDBN
	static int L[16] = {
	    5,0,1,5, 2,5,5,5, 3,5,5,5, 5,5,5,5
	};

	// convert from sam base to acgt*n order.
	base = L[base];
	if (p[n].is_del) base = 4;

	double MM, __, _M, qe;

	// Correction for mapping quality.  Maybe speed up via lookups?
	// Cannot nullify mapping quality completely.  Lots of (true)
	// SNPs means low mapping quality.  (Ideally need to know
	// hamming distance to next best location.)

	if (flags & CONS_MQUAL) {
	    double _p = mqual_pow[qual];
	    double _m = mqual_pow[b->core.qual];

	    //printf("%c %d -> %d, %f %f\n", "ACGT*N"[base], qual, (int)(-TENOVERLOG10 * log(1-(_m * _p + (1 - _m)/4))), _p, _m);
	    qual = ph_log(1-(_m * _p + (1 - _m)/4));
	}

	/* Quality 0 should never be permitted as it breaks the math */
	if (qual < 1)
	    qual = 1;

	__ = p__[stech][qual];
	MM = pMM[stech][qual] - __;
	_M = p_M[stech][qual] - __;

	if (flags & CONS_DISCREP) {
	    qe = q2p[qual];
	    sumsE += qe;
	    sumsC[base] += 1 - qe;
	}

	counts[base]++;

	switch (base) {
	case 0:
	    S[0] += MM; S[1 ]+= _M; S[2 ]+= _M; S[3 ]+= _M; S[4 ]+= _M;
	    break;

	case 1:
	    S[1 ]+= _M; S[5 ]+= MM; S[6 ]+= _M; S[7 ]+= _M; S[8 ]+= _M;
	    break;

	case 2:
	    S[2 ]+= _M; S[6 ]+= _M; S[9 ]+= MM; S[10]+= _M; S[11]+= _M; 
	    break;

	case 3:
	    S[3 ]+= _M; S[7 ]+= _M; S[10]+= _M; S[12]+= MM; S[13]+= _M; 
	    break;

	case 4:
	    S[4 ]+= _M; S[8 ]+= _M; S[11]+= _M; S[13]+= _M; S[14]+= MM;
	    break;

	case 5: /* N => equal weight to all A,C,G,T but not a pad */
	    S[0] += MM; S[1 ]+= MM; S[2 ]+= MM; S[3 ]+= MM; S[4 ]+= _M;
	                S[5 ]+= MM; S[6 ]+= MM; S[7 ]+= MM; S[8 ]+= _M;
			            S[9 ]+= MM; S[10]+= MM; S[11]+= _M; 
				                S[12]+= MM; S[13]+= _M; 
	    break;
	}

	depth++;
    }


    /* and speculate */
    {
	double shift, max, max_het, norm[15];
	int call = 0, het_call = 0, ph;
	double tot1, tot2;

	/*
	 * Scale numbers so the maximum score is 0. This shift is essentially 
	 * a multiplication in non-log scale to both numerator and denominator,
	 * so it cancels out. We do this to avoid calling exp(-large_num) and
	 * ending up with norm == 0 and hence a 0/0 error.
	 *
	 * Can also generate the base-call here too.
	 */
	shift = -DBL_MAX;
	max = -DBL_MAX;
	max_het = -DBL_MAX;

	for (j = 0; j < 15; j++) {
	    S[j] += lprior15[j];
	    if (shift < S[j])
		shift = S[j];

	    /* Only call pure AA, CC, GG, TT, ** for now */
	    if (j != 0 && j != 5 && j != 9 && j != 12 && j != 14) {
		if (max_het < S[j]) {
		    max_het = S[j];
		    het_call = j;
		}
		continue;
	    }

    if (max < S[j]) {
		max = S[j];
		call = j;
	    }
	}

	/*
	 * Shift and normalise.
	 * If call is, say, b we want p = b/(a+b+c+...+n), but then we do
	 * p/(1-p) later on and this has exceptions when p is very close
	 * to 1.
	 *
	 * Hence we compute b/(a+b+c+...+n - b) and
	 * rearrange (p/norm) / (1 - (p/norm)) to be p/norm2.
	 */
	for (j = 0; j < 15; j++) {
	    S[j] -= shift;
	    double e = fast_exp(S[j]);
	    S[j] = (S[j] > min_e_exp) ? e : DBL_MIN;
	    norm[j] = 0;
	}

	tot1 = tot2 = 0;
	for (j = 0; j < 15; j++) {
	    norm[j]    += tot1;
	    norm[14-j] += tot2;
	    tot1 += S[j];
	    tot2 += S[14-j];
	}

	/* And store result */
	if (depth && depth != counts[5] /* all N */) {
	    double m;

	    cons->depth = depth;

	    cons->call     = map_sing[call];
	    if (norm[call] == 0) norm[call] = DBL_MIN;
	    ph = ph_log(norm[call]) + .5;
	    cons->phred = ph > 255 ? 255 : (ph < 0 ? 0 : ph);
	    //cons->call_prob1 = norm[call]; // p = 1 - call_prob1

	    cons->het_call = map_het[het_call];
	    if (norm[het_call] == 0) norm[het_call] = DBL_MIN;
	    ph = TENLOG2OVERLOG10 * (fast_log2(S[het_call]) - fast_log2(norm[het_call])) + .5;

	    cons->het_phred = ph;
	    //cons->het_prob_n = S[het_call]; // p = prob_n / prob_d
	    //cons->het_prob_d = norm[het_call];

	    /* Compute discrepancy score */
	    if (flags & CONS_DISCREP) {
		m = sumsC[0]+sumsC[1]+sumsC[2]+sumsC[3]+sumsC[4];
		double c;
		if (cons->het_phred > 0)
		    c = sumsC[cons->het_call%5] + sumsC[cons->het_call/5];
		else
		    c = sumsC[cons->call];;
		cons->discrep = (m-c)/sqrt(m);
//		printf("Discrep = %f,  %f %f %f %f %f\n", cons->discrep,
//		       sumsC[0], sumsC[1], sumsC[2], sumsC[3], sumsC[4]);
//		if (cons->discrep > 1)
//		    printf("XYZZY\n");
	    }
	} else {
	    cons->call = 5; /* N */
	    cons->het_call = 0;
	    cons->het_phred = 0;
	    cons->phred = 0;
	    cons->depth = 0;
	    cons->discrep = 0;
	}
    }

    return 0;
}

// Applies a basic P-block algorithm to a string of quality values.
// All blocks of qualities within +/- level get replaced by a representative
// value such that the delta is within 'level'.
// The original paper had a more complex method, but this is just (min+max)/2.
void pblock(bam1_t *b, int level) {
    if (!b)
	return;

    int len = b->core.l_qseq, i, j, qmin = INT_MAX, qmax = INT_MIN;
    int last_qmin = 0, last_qmax = 0, mid;
    uint8_t *qual = bam_get_qual(b);

    level *= 2;

    for (i = j = 0; i < len; i++) {
	if (qmin > qual[i])
	    qmin = qual[i];
	if (qmax < qual[i])
	    qmax = qual[i];
	if (qmax - qmin > level) {
	    mid = (last_qmin + last_qmax) / 2;
	    memset(qual+j, mid, i-j);
	    qmin = qmax = qual[i];
	    j = i;
	}
	last_qmin = qmin;
	last_qmax = qmax;
    }

    mid = (last_qmin + last_qmax) / 2;
    memset(qual+j, mid, i-j);
}

// // Quantise qualities by bin[] array.  Alternative to pblock above
// void qblock(bam1_t *b) {
//     if (!b)
// 	return;
// 
//     int len = b->core.l_qseq, i;
//     uint8_t *qual = bam_get_qual(b);
// 
//     for (i = 0; i < len; i++)
// 	qual[i] = bin[qual[i]];
// }

//-----------------------------------------------------------------------------
// Tree of bam objects, sorted by chromosome & position

typedef struct bam_sorted_item {
    RB_ENTRY(bam_sorted_item) link;
    bam1_t *b;
    uint64_t id;
    int end_pos;
    int keep_qual;
} bam_sorted_item;

RB_HEAD(bam_sort, bam_sorted_item);

typedef struct bam_sort bam_sorted_list;

static int bam_item_cmp(bam_sorted_item *b1, bam_sorted_item *b2) {
    int d;
    if (b2->b->core.tid == -1)
	return -1;
    if ((d = b1->b->core.tid - b2->b->core.tid))
	return d;
    if ((d = b1->b->core.pos - b2->b->core.pos))
	return d;
    return b1->id - b2->id;
}

RB_PROTOTYPE(bam_sort, bam_sorted_item, link, bam_item_cmp);
RB_GENERATE(bam_sort, bam_sorted_item, link, bam_item_cmp);

bam_sorted_list *bam_sorted_list_new(void) {
    bam_sorted_list *bl = calloc(1, sizeof(bam_sorted_list));
    if (!bl)
        return NULL;
        
    RB_INIT(bl);

    return bl;
}

void bam_sorted_list_destroy(bam_sorted_list *bl) {
    bam_sorted_item *node, *next;
    
    if (!bl)
        return;

    for (node = RB_MIN(bam_sort, bl); node; node = next) {
	next = RB_NEXT(bam_sort, bl, node);
	RB_REMOVE(bam_sort, bl, node);
	free(node);
    }

    free(bl);
}

bam_sorted_item *bam_sorted_item_new(void) {
    return calloc(1, sizeof(bam_sorted_item));
}

/*
 * Inserts a bam structure into the bam_sorted_list, in positional
 * order. 
 *
 * Returns the newly created bam_sorted_list element on success.
 *         NULL on failure.
 */
bam_sorted_item *insert_bam_list2(bam_sorted_list *bl, bam_sorted_item *ele) {
    RB_INSERT(bam_sort, bl, ele);

    return ele;
}

bam_sorted_item *insert_bam_list_id(bam_sorted_list *bl, bam1_t *b, uint64_t id) {
    bam_sorted_item *ele = bam_sorted_item_new();
    if (!ele)
	return NULL;
    
    ele->b = b;
    ele->id = id;
    ele->end_pos = bam_endpos(ele->b);

    return insert_bam_list2(bl, ele);
}

static uint64_t global_id = 0;
bam_sorted_item *insert_bam_list(bam_sorted_list *bl, bam1_t *b) {
    return insert_bam_list_id(bl, b, global_id++);
}

/*
 * Removes an item from the bam_sorted_list.
 */
void remove_bam_list(bam_sorted_list *bl, bam_sorted_item *ele) {
    RB_REMOVE(bam_sort, bl, ele);

    free(ele);
}

//-----------------------------------------------------------------------------
// Copied from htslib/sam.c.
// TODO: we need a proper interface to find the length of an aux tag,
// or at the very make exportable versions of these in htslib.
static inline int aux_type2size(uint8_t type)
{
    switch (type) {
    case 'A': case 'c': case 'C':
        return 1;
    case 's': case 'S':
        return 2;
    case 'i': case 'I': case 'f':
        return 4;
    case 'd':
        return 8;
    case 'Z': case 'H': case 'B':
        return type;
    default:
        return 0;
    }
}

// Copied from htslib/sam.c.
static inline uint8_t *skip_aux(uint8_t *s)
{
    int size = aux_type2size(*s); ++s; // skip type
    uint32_t n;
    switch (size) {
    case 'Z':
    case 'H':
        while (*s) ++s;
        return s + 1;
    case 'B':
        size = aux_type2size(*s); ++s;
        memcpy(&n, s, 4); s += 4;
        return s + size * n;
    case 0:
        abort();
        break;
    default:
        return s + size;
    }
}

void purge_tags(cram_lossy_params *settings, bam1_t *b) {
    if (settings->aux_whitelist) {
        uint8_t *s_from, *s_to;
        auxhash_t h = settings->aux_whitelist;

        s_from = s_to = bam_get_aux(b);
        while (s_from < b->data + b->l_data) {
            int x = (int)s_from[0]<<8 | s_from[1];
            uint8_t *s = skip_aux(s_from+2);

            if (kh_get(aux_exists, h, x) != kh_end(h) ) {
                if (s_to != s_from) memmove(s_to, s_from, s - s_from);
                s_to += s - s_from;
            }
            s_from = s;
        }
        b->l_data = s_to - b->data;

    } else if (settings->aux_blacklist) {
        uint8_t *s_from, *s_to;
        auxhash_t h = settings->aux_blacklist;

        s_from = s_to = bam_get_aux(b);
        while (s_from < b->data + b->l_data) {
            int x = (int)s_from[0]<<8 | s_from[1];
            uint8_t *s = skip_aux(s_from+2);

            if (kh_get(aux_exists, h, x) == kh_end(h) ) {
                if (s_to != s_from) memmove(s_to, s_from, s - s_from);
                s_to += s - s_from;
            }
            s_from = s;
        }
        b->l_data = s_to - b->data;
    }

    if (settings->BD_low || settings->BD_mid || settings->BD_high) {
        uint8_t *t = bam_get_aux(b);
        while (t < b->data + b->l_data) {
	    if (t[0] == 'B' && t[1] == 'D') {
		uint8_t *c = t+2;
		while(*++c) {
		    *c = (*c >= settings->BD_mid)
			? settings->BD_high
			: settings->BD_low;
		}
	    }
	    t = skip_aux(t+2);
        }
    }

    if (settings->BI_low || settings->BI_mid || settings->BI_high) {
        uint8_t *t = bam_get_aux(b);
        while (t < b->data + b->l_data) {
	    if (t[0] == 'B' && t[1] == 'I') {
		uint8_t *c = t+2;
		while(*++c) {
		    *c = (*c >= settings->BI_mid)
			? settings->BI_high
			: settings->BI_low;
		}
	    }
	    t = skip_aux(t+2);
        }
    }
}

//-----------------------------------------------------------------------------
typedef struct {
    // sam, bam AND hts! All the weasles in a single bag :)
    samFile *fp;
    bam_hdr_t *header;
    hts_itr_t *iter;
    bam_sorted_list *bl, *b_hist;
    bam1_t *b_unmap;
    int64_t count_in, count_out;
} pileup_cd;

int flush_bam_list(pileup_cd *cd, cram_lossy_params *p, bam_sorted_list *bl,
		   int before_tid, int before, samFile *out, bam_hdr_t *header) {
    bam_sorted_item *bi, *next;

    // Sanity check
    int last_pos = 0;
    RB_FOREACH(bi, bam_sort, bl) {
	if (!(bi->b->core.flag & BAM_FUNMAP)) {
	    assert(bi->b->core.pos >= last_pos);
	    last_pos = bi->b->core.pos;
	}
    }

    for (bi = RB_MIN(bam_sort, bl); bi; bi = next) {
	next = RB_NEXT(bam_sort, bl, bi);
	if (bi->end_pos >= before ||
	    (bi->b->core.tid >= 0 && bi->b->core.tid >= before_tid))
	    break;

	purge_tags(p, bi->b);

	// Correct qualities
	int x;
	uint8_t *qual = bam_get_qual(bi->b);
	for (x = 0; x < bi->b->core.l_qseq; x++) {
	    if (qual[x] & 0x80)
		qual[x] &= ~0x80;
	}
	cd->count_out++;
	if (p->pblock)
	    pblock(bi->b, p->pblock);
	if (sam_write1(out, header, bi->b) < 0)
	    return -1;
	bam_destroy1(bi->b);
	remove_bam_list(bl, bi);
    }

    return 0;
}

//-----------------------------------------------------------------------------
// Main pileup iterator

int pileup_callback(void *vp, bam1_t *b) {
    pileup_cd *cd = (pileup_cd *)vp;
    int ret = (cd->iter)
	? sam_itr_next(cd->fp, cd->iter, b)
	: sam_read1(cd->fp, cd->header, b);

    if (ret >= 0) {
	cd->count_in++;

	// Unmapped chromosome => end of pileup.  Record the read we've
	// already read and then feign EOF so we can handle these outside
	// of the pileup interface.
	if (b->core.tid == -1) {
	    cd->b_unmap = bam_dup1(b);
	    return -1;
	}

	int unmap = b->core.flag & BAM_FUNMAP;
	if (!unmap) {
	    // Mapped reads that have no matching ref location, eg
	    // they are entirely insertion, will not appear in any pileup
	    // column.  Therefore treat them as unmapped to avoid errors.
	    uint32_t *cig = bam_get_cigar(b);
	    int i, n = b->core.n_cigar;
	    for (i = 0; i < n; i++)
		if (bam_cigar_type(bam_cigar_op(cig[i])) & 2)
		    break;
	    if (i == n) {
#ifdef DEBUG
		printf("Note: %s mapped at #%d,%d but has no cigar ref op!\n",
		       bam_get_qname(b), b->core.tid, b->core.pos);
#endif
		unmap = 1;
	    }
	}

	insert_bam_list(unmap ? cd->b_hist : cd->bl, bam_dup1(b));
    }

    return ret;
}

// Turns an absolute reference position into a relative query position within the seq.
int ref2query_pos(bam1_t *b, int pos) {
     uint32_t *cig = bam_get_cigar(b);
     int i, n = b->core.n_cigar;
     int p = b->core.pos, q = 0;

     for (i = 0; i < n; i++) {
 	int op = bam_cigar_op(cig[i]);
 	int op_len = bam_cigar_oplen(cig[i]);
	if (p + ((bam_cigar_type(op) & 2) ? op_len : 0) < pos) {
	    if (bam_cigar_type(op) & 1) // query
		q += op_len;
	    if (bam_cigar_type(op) & 2) // ref
		p += op_len;
	    continue;
	}

	if (bam_cigar_type(op) & 1) // query
	    q += (pos - p); // consume partial op_len

	return q >= 0 ? q : 0;
     }

     return q;
}


// // Masks any base within INDEL_ADD bases of a soft-clip. Hard too?
// void mask_clips(bam1_t *b) {
//     uint32_t *cig = bam_get_cigar(b);
//     int i, n = b->core.n_cigar, q;
// 
//     for (i = q = 0; i < n; i++) {
// 	int op = bam_cigar_op(cig[i]);
// 	int oplen = bam_cigar_oplen(cig[i]);
// 
// 	if (op == BAM_CSOFT_CLIP) {
// 	    uint8_t *qual = bam_get_qual(b);
// 	    int x, x_s = q-INDEL_ADD, x_e = q+oplen+INDEL_ADD;
// 	    if (x_s < 0) x_s = 0;
// 	    if (x_e > b->core.l_qseq) x_e = b->core.l_qseq;
// 	    for (x = x_s; x < x_e; x++)
// 		qual[x] |= 0x80;
// 	}
// 	if (bam_cigar_type(op) & 1)
// 	    q += oplen;
//     }
// }

// Converts a position on a query sequence to a position on the reference.
int bam_qpos2rpos(bam1_t *b, int qpos) {
    const uint32_t *cigar = bam_get_cigar(b);
    int k, rpos = b->core.pos, aqpos = 0; // accumulated qpos
    for (k = 0; k < b->core.n_cigar && aqpos < qpos; k++) {
	if (bam_cigar_type(bam_cigar_op(cigar[k]))&2) {
	    if (bam_cigar_oplen(cigar[k]) <= qpos - aqpos)
		rpos += bam_cigar_oplen(cigar[k]);
	    else
		rpos += qpos - aqpos;
	}
	if (bam_cigar_type(bam_cigar_op(cigar[k]))&1)
	    aqpos += bam_cigar_oplen(cigar[k]);
    }
    return rpos;
}

// Extend min/max reference positions based on any STRs at apos/rpos (abs/rel)
// into seq b.  This is used to find the extents over which we may wish to
// preserve scores given something important is going on at apos (typically
// indel).
int mask_LC_regions(cram_lossy_params *p, int is_indel, bam1_t *b,
		    int apos, int rpos, int *min_pos, int *max_pos) {
    int len = b->core.l_qseq;
    char *seq = malloc(len);
    int i;

    if (!seq)
	return -1;

    for (i = 0; i < len; i++)
	seq[i] = seq_nt16_str[bam_seqi(bam_get_seq(b), i)];

    rep_ele *reps = find_STR(seq, len, 0), *elt, *tmp;

    DL_FOREACH_SAFE(reps, elt, tmp) {
	if (is_indel && 
	    !(rpos+p->iSTR_add >= elt->start && rpos-p->iSTR_add <= elt->end)) {
	    //fprintf(stderr, "SKIP rpos %d, apos %d:\t%2d .. %2d %.*s\n",
	    //	    rpos, apos,
	    //	    elt->start, elt->end,
	    //	    elt->end - elt->start+1, &seq[elt->start]);
	    DL_DELETE(reps, elt);
	    free(elt);
	    continue;
	} else if (!is_indel &&
		   !(rpos+p->sSTR_add >= elt->start && rpos-p->sSTR_add <= elt->end)) {
	    DL_DELETE(reps, elt);
	    free(elt);
	    continue;
	}
	
	//fprintf(stderr, "%d:\t%2d .. %2d %.*s\n",
	//	apos,
	//	elt->start, elt->end,
	//	elt->end - elt->start+1, &seq[elt->start]);

	elt->start = bam_qpos2rpos(b, elt->start);
	elt->end   = bam_qpos2rpos(b, elt->end);
	if (is_indel) {
	    if (*min_pos > elt->start)
		*min_pos = elt->start;
	    if (*max_pos < elt->end)
		*max_pos = elt->end;
	} else {
	    if (*min_pos > elt->start)
		*min_pos = elt->start;
	    if (*max_pos < elt->end)
		*max_pos = elt->end;
	}

	DL_DELETE(reps, elt);
	free(elt);
    }    

    free(seq);
    return 0;
}

static int count_het_qual_A = 0;
static int count_het_qual_B = 0;
static int count_hom_qual_A = 0;
static int count_hom_qual_B = 0;
static int count_het_A = 0;
static int count_het_B = 0;
static int count_hom_A = 0;
static int count_hom_B = 0;
static int count_discrep_A = 0;
static int count_discrep_B = 0;
static int count_diff = 0;
static int count_indel = 0;
static int count_indel_qual = 0;

static int64_t count_columns = 0;
static int64_t count_low_mqual_perc = 0;
static int64_t count_clip_perc = 0;
static int64_t count_ins_len_perc = 0;
static int64_t count_indel_ov_perc = 0;
static int64_t count_over_depth = 0;

static char *tid_name(bam_hdr_t *h, int tid) {
    return h->target_name[tid];
}

// FIXME: this needs breaking down; it's become far too bloated
int transcode(cram_lossy_params *p, samFile *in, samFile *out,
	      bam_hdr_t *header, hts_itr_t *h_iter) {
    bam_plp_t p_iter;
    int tid, pos, last_tid = -2;
    int n_plp;
    const bam_pileup1_t *plp;
    pileup_cd cd = {0};
    bam_sorted_list *bl = bam_sorted_list_new();
    bam_sorted_list *b_hist = bam_sorted_list_new();
    int str_snp = (p->sSTR_add || p->sSTR_mul);
    int counter = 0;

    cd.fp = in;
    cd.header = header;
    cd.iter = h_iter;
    cd.bl = bl;
    cd.b_hist = b_hist;
    cd.count_in = 0;
    cd.count_out = 0;

    p_iter = bam_plp_init(pileup_callback, &cd);
    bam_plp_set_maxcnt(p_iter, INT_MAX);
    int min_pos = INT_MAX, max_pos = 0;
    int min_pos2 = INT_MAX, max_pos2 = 0;

    int64_t total_depth = 0, total_col = 0;

    while ((plp = bam_plp_auto(p_iter, &tid, &pos, &n_plp))) {
	int i, preserve = 0, indel = 0;
	unsigned char base;
	int left_most = n_plp ? plp[0].b->core.pos : 0;

	count_columns++;

	if (tid != last_tid) {
	    // Ensure b_hist is only per chromosome
	    if (flush_bam_list(&cd, p, b_hist, tid, INT_MAX, out, header) < 0)
		return -1;
	    last_tid = tid;
	    min_pos = INT_MAX, max_pos = 0;
	    min_pos2 = INT_MAX, max_pos2 = 0;
	    total_depth = 0;
	    total_col = 0;
	}

	total_depth += n_plp;
	total_col++;

	if (n_plp > MAX_DEPTH) {
	    if (p->verbose > 1)
		fprintf(stderr, "Excessive depth at tid %d, pos %d, depth %d\n", tid, pos, n_plp);
	    if (p->bed_fp)
		fprintf(p->bed_fp, "%s\t%d\t%d\tVDEEP\n",
			tid_name(header, tid), MAX(pos-BED_DIST,0), pos+BED_DIST);
	    goto too_deep;
	}

	if (counter++ == 100000) {
	    if (p->verbose)
		fprintf(stderr, "Processing %s:%d\n", tid_name(header,tid), pos);
	    counter = 0;
	}

	if (pos > max_pos2) {
	    min_pos2 = min_pos = INT_MAX;
	    max_pos2 = max_pos = 0; // beyond region to preserve quals.
	}

	if (h_iter) {
	    if (pos < h_iter->beg)
		continue;
	    if (pos >= h_iter->end)
		break;
	}

	consensus_t cons_g5_A, cons_g5_B;
	int call1 = 0, call2 = 0; // samtools numerical =ACMGRSVTWYHKDBN

	if (p->min_qual_A != 0) {
	    calculate_consensus_pileup(CONS_DISCREP,
				       plp, n_plp, &cons_g5_A);
	    if (cons_g5_A.het_phred > 0) {
		call1 = 1<<(cons_g5_A.het_call / 5);
		call2 = 1<<(cons_g5_A.het_call % 5);
	    } else {
		call1 = call2 = 1<<cons_g5_A.call;
	    }
	}

	if (p->min_qual_B != 0) {
	    calculate_consensus_pileup(CONS_DISCREP | CONS_MQUAL,
				       plp, n_plp, &cons_g5_B);
	    if (cons_g5_B.het_phred > 0) {
		call1 = 1<<(cons_g5_B.het_call / 5);
		call2 = 1<<(cons_g5_B.het_call % 5);
	    } else {
		call1 = call2 = 1<<cons_g5_B.call;
	    }
	}

#ifdef DEBUG
	printf("Depth %d\t%d\t%d\t", tid, pos+1, n_plp);
#endif

#ifdef DEBUG
	if (p->min_qual_A != 0) {
	    if (cons_g5_A.het_phred > 0) {
		printf("%c/%c %4d\t",
		       "ACGT*"[cons_g5_A.het_call / 5], "ACGT*"[cons_g5_A.het_call % 5],
		       (int)cons_g5_A.het_phred);
	    } else {
		printf("%c   %4d\t",
		       "ACGT*N"[cons_g5_A.call],
		       cons_g5_A.phred);
	    }
	}

	if (p->min_qual_B != 0) {
	    if (cons_g5_B.het_phred > 0) {
		printf("%c/%c %4d\t",
		       "ACGT*"[cons_g5_B.het_call / 5], "ACGT*"[cons_g5_B.het_call % 5],
		       (int)cons_g5_B.het_phred);
	    } else {
		printf("%c   %4d\t",
		       "ACGT*N"[cons_g5_B.call],
		       cons_g5_B.phred);
	    }
	}
#endif

	int hA = 0, sA = 0, hB = 0, sB = 0;
	if (p->min_qual_A) {
	    hA = cons_g5_A.het_phred > 0
		? cons_g5_A.het_call
		: cons_g5_A.call * 5 + cons_g5_A.call;
	    sA = cons_g5_A.het_phred > 0
		? cons_g5_A.het_phred
		: cons_g5_A.phred;
	}
	if (p->min_qual_B) {
	    hB = cons_g5_B.het_phred > 0
		? cons_g5_B.het_call
		: cons_g5_B.call * 5 + cons_g5_B.call;
	    sB = cons_g5_B.het_phred > 0
		? cons_g5_B.het_phred
		: cons_g5_B.phred;
	}

	if (p->min_qual_A && p->min_qual_B && hA != hB) count_diff++;
	if (p->min_qual_A) {
	    if (cons_g5_A.het_phred > 0) {
		count_het_A++;
		if (sA < p->min_qual_A)
		    count_het_qual_A++;
	    } else {
		count_hom_A++;
		if (sA < p->min_qual_A)
		    count_hom_qual_A++;
	    }
	    if (cons_g5_A.discrep >= p->min_discrep_A)
		count_discrep_A++;
	}
	if (p->min_qual_B) {
	    if (cons_g5_B.het_phred > 0) {
		count_het_B++;
		if (sB < p->min_qual_B)
		    count_het_qual_B++;
	    } else {
		count_hom_B++;
		if (sB < p->min_qual_B)
		    count_hom_qual_B++;
	    }
	    if (cons_g5_B.discrep >= p->min_discrep_B)
		count_discrep_B++;
	}

	preserve = ((p->min_qual_A && p->min_qual_B && hA != hB) ||
		    (p->min_qual_A && sA < p->min_qual_A) ||
		    (p->min_qual_B && sB < p->min_qual_B));
	preserve |= ((p->min_qual_A && cons_g5_A.discrep >= p->min_discrep_A) ||
		     (p->min_qual_B && cons_g5_B.discrep >= p->min_discrep_B));

#ifdef DEBUG
	if (preserve) {
	    printf("*\t");
	} else {
	    printf("\t");
	}
#endif

	// Check average mapping quality.
	int had_indel = 0, had_indel_Q = 0;
	int low_mq_count = 0, keep_qual = 0;
	for (i = 0; i < n_plp; i++) {
	    low_mq_count += (plp[i].b->core.qual <= p->min_mqual);
	    if (plp[i].indel || plp[i].is_del)
		had_indel = 1;
	}

	keep_qual = low_mq_count > p->low_mqual_perc * (n_plp + .01);
	count_low_mqual_perc += keep_qual;

	// Check for unexpectedly deep regions.
	//printf("%d: %lld %lld %d\n", pos, total_depth, total_col, total_depth/total_col);
	if (n_plp*(total_col+1) > p->over_depth * (total_depth+1)) {
	    //fprintf(stderr, "%d %d\tUnexpectedly high depth: %d vs %d\n",
	    //	    tid, pos, n_plp, (int)(total_depth / (total_col+1)));
	    if (p->bed_fp)
		fprintf(p->bed_fp, "%s\t%d\t%d\tDEEP\n",
			tid_name(header, tid), MAX(pos-BED_DIST,0), pos+BED_DIST);
	    keep_qual = 1;
	    count_over_depth++;
	}

	// Keep average depth to within the last few Mb only.
	if (total_col > 1024*1024) {
	    total_col   >>= 1;
	    total_depth >>= 1;
	}


	// Look for indel and compute STR extents.
	int indel_sz = 0;
	int indel_depth[101];
	indel_depth[0] = 0;
	int clipped = 0, n_overlap = 0;
	for (i = 0; i < n_plp; i++) {
	    int is_indel = (plp[i].indel || plp[i].is_del);

	    if ((plp[i].is_head && plp[i].qpos > 0) ||
		(plp[i].is_tail && plp[i].qpos+1 < plp[i].b->core.l_qseq))
		clipped++;

	    if (!plp[i].is_tail && !plp[i].is_head)
		n_overlap++;

	    if (!plp[i].is_head && !plp[i].is_tail && (plp[i].indel > 0 || had_indel)) {
		while (indel_sz < plp[i].indel && indel_sz < 100)
		    indel_depth[++indel_sz] = 0;
		if (plp[i].indel >= 0)
		    indel_depth[MIN(plp[i].indel, 99)]++;
	    }

	    // For indels, expand to cover low-complexity regions.
	    // (Also for SNPs if required via -i option.)
	    // -q0 or -Q0 (min_qual_A/B) disables that cons algorithm.
	    if ((is_indel || (str_snp && preserve))
		&& ((p->min_qual_A && sA < p->min_indel_A) ||
		    (p->min_qual_B && sB < p->min_indel_B))) {

		if (is_indel) had_indel_Q++;

		if (is_indel) {
		    if (indel < ABS(plp[i].indel) + plp[i].is_del)
			indel = ABS(plp[i].indel) + plp[i].is_del;
		} else {
		    indel = 1;
		}

		if (mask_LC_regions(p, is_indel, plp[i].b, pos,
				    plp[i].qpos+1, &min_pos, &max_pos) < 0)
		    return -1;
		if (mask_LC_regions(p, is_indel, plp[i].b, pos+indel,
				    plp[i].qpos+1, &min_pos, &max_pos) < 0)
		    return -1;

		if (min_pos > pos) min_pos = pos;
		if (max_pos < pos) max_pos = pos;

		// Extra growth, paranoia.
		if (is_indel) {
		    min_pos2 = MIN(min_pos2,
				   pos - (pos-min_pos)*p->iSTR_mul - p->iSTR_add);
		    max_pos2 = MAX(max_pos2,
				   pos + (max_pos-pos)*p->iSTR_mul + p->iSTR_add);
		} else {
		    min_pos2 = MIN(min_pos2,
				   pos - (pos-min_pos)*p->sSTR_mul - p->sSTR_add);
		    max_pos2 = MAX(max_pos2,
				   pos + (max_pos-pos)*p->sSTR_mul + p->sSTR_add);
		}
	    }

	    //if (min_pos != INT_MAX || indel)
	    //	fprintf(stderr, "%d..%d: Mask region = %d..%d\n", pos-1, pos+indel+1, min_pos, max_pos);
	}
	if (had_indel) count_indel++;
	if (had_indel_Q) count_indel_qual++;

	if ((clipped - 1.0) >= p->clip_perc * n_overlap) {
	    if (p->verbose > 1)
		fprintf(stderr, "%s %d\tUnexpected high clip rate, %d of %d\n",
			tid_name(header,tid), pos, clipped, n_overlap);
	    if (p->bed_fp)
		fprintf(p->bed_fp, "%s\t%d\t%d\tCLIP\n",
			tid_name(header, tid), MAX(pos-BED_DIST,0), pos+BED_DIST);
	    keep_qual = 1;
	    count_clip_perc++;
	}

	// Over the span of an indel, the ratio of top to total should
	// be consistently 1 or 2 things, similarly the depth.
	if (indel_sz) {
	    //printf("Indel at %d, max size %d\n", pos, indel_sz);
	    //int qv1 = 0, qv2 = 0;
	    int qd1 = 0, qd2 = 0;
	    int indel_overlap = 0;
	    for (i = 0; i <= indel_sz && i < 100; i++) {
		if (!indel_depth[i])
		    continue;
	    	//printf("%3d\t%2d\n", i, indel_depth[i]);
		indel_overlap += indel_depth[i];
		if (qd1 < indel_depth[i]) {
		    qd2 = qd1; //qv2 = qv1;
		    qd1 = indel_depth[i];
		    //qv1 = i;
		} else if (qd2 < indel_depth[i]) {
		    qd2 = indel_depth[i];
		    //qv2 = i;
		}
	    }
	    //printf("Top 2 = %d x %d,  %d x %d, out of %d, ov/n_plp=%f\n", qv1, qd1, qv2, qd2, indel_overlap, (double)indel_overlap / n_plp);
	    
	    if ((indel_overlap - qd1 - qd2) > p->ins_len_perc * (indel_overlap + .1)) {
		if (p->verbose > 1)
		    fprintf(stderr, "%s %d\tSuspect indel, depth %d / %d, common %d+%d\n",
			    tid_name(header,tid), pos, n_plp, indel_overlap, qd1, qd2);
		if (p->bed_fp)
		    fprintf(p->bed_fp, "%s\t%d\t%d\tINDEL_LEN\n",
			    tid_name(header, tid), MAX(pos-BED_DIST,0), pos+BED_DIST);
		keep_qual = 1;
		count_ins_len_perc++;
	    }

	    if ((double)indel_overlap < p->indel_ov_perc * n_plp) {
		if (p->bed_fp)
		    fprintf(p->bed_fp, "%s\t%d\t%d\tINDEL_COVERAGE\n",
			    tid_name(header, tid), MAX(pos-BED_DIST,0), pos+BED_DIST);
		if (p->verbose > 1)
		    fprintf(stderr, "%s %d\tSuspect drop in indel overlap %d vs %d\n",
			    tid_name(header,tid), pos, indel_overlap, n_plp);
		keep_qual = 1;
		count_indel_ov_perc++;
	    }
	}


	// Mask individual bases
	bam_sorted_item *bi = RB_MIN(bam_sort, bl);
	for (i = 0; i < n_plp; i++) {
	    // b is unedited bam in pileup.
	    // b2 is edited bam in bam_sorted_list, for output
	    bam1_t *b = plp[i].b, *b2 = NULL;
	    uint8_t *qual;

	    assert(bi);

	    b2 = bi->b;

	    // We have an assumption that all mapped reads appear in the pileup
	    // somewhere.
	    //
	    // The caveat is any read that doesn't have a M, =, X or D cigar
	    // op (eg "5S20I").  We treat these as unmapped in pileup_callback
	    // to avoid this assertion failing.

	    // Assert b2 & b are same object.
	    assert(b2 && strcmp(bam_get_qname(b), bam_get_qname(b2)) == 0);

	    // Any column in this alignment has the chance to set keep_qual.
	    // Setting it means all qualities for the overlapping sequences
	    //are retained.
	    if (keep_qual) bi->keep_qual = 1;

	    bi = RB_NEXT(bam_sort, bl, bi);

	    // First time for this seq, handle low mqual if desired.
	    if (plp[i].is_head) {
		//mask_clips(b2);
		if (b->core.qual <= p->min_mqual) {
		    int x;
		    for (x = 0; x < b->core.l_qseq; x++)
			bam_get_qual(b2)[x] |= 0x80;
		}
	    }

	    qual = &bam_get_qual(b2)[plp[i].qpos];
	    base = bam_seqi(bam_get_seq(b), plp[i].qpos);
#ifdef DEBUG
	    putchar(seq_nt16_str[base]);
#endif

	    if (indel) {
		// New indel, so backfill to low-complexity range.
		int x;
		//fprintf(stderr, "pos %d indel %d, min/max_pos %d/%d, core pos %d, qpos %d => %d\n",
		//	pos, indel, min_pos2, max_pos2, b->core.pos, plp[i].qpos, ref2query_pos(b, min_pos2));
		for (x = ref2query_pos(b, min_pos2); x <= plp[i].qpos; x++) {
		    //bam_get_qual(b2)[x] = 11 | 0x80;
		    bam_get_qual(b2)[x] = bam_get_qual(b)[x] | 0x80;
		}
	    }
	    if (min_pos != INT_MAX) { // or max_pos != 0; ie not reset yet
		// not indel at this pos, but still in range.
		*qual = bam_get_qual(b)[plp[i].qpos] | 0x80;
	    }

	    if (preserve)
		*qual |= 0x80;

	    if (!keep_qual && p->softclip) {
		if (plp[i].is_head) {
		    int x;
		    for (x = 1; x <= plp[i].qpos; x++)
			qual[-x] = bin2[qual[-x]];
		} else if (plp[i].is_tail) {
		    int x;
		    for (x = b->core.l_qseq - plp[i].qpos -1; x>0; x--)
			qual[x] = bin2[qual[x]];
		}
	    }

	    if (!(*qual & 0x80)) {
		// no need to preserve quals
		if (base == call1 || base == call2) {
		    // base matches hom or het call; increase confidence to fixed amount
		    *qual = p->qhigh;
		} else if (p->reduce_qual) {
		    // base mismatches call, but call is confident. => reduce qual.
		    // bases in the consensus that don't match the expected heterozygous call.
		    if (p->binary_qual)
			*qual = bin2[*qual];
		    else
			*qual = p->qlow;
		}
	    }
	}

#ifdef DEBUG
	printf("\n");
#endif

    too_deep:

	// Push any finished sequence to the b_hist list,
	// clearing the qual mask bit as we go.
	bi = RB_MIN(bam_sort, bl);
	for (i = 0; i < n_plp; i++) {
	    bam1_t *b2 = bi->b;

	    if (!plp[i].is_tail) {
		bi = RB_NEXT(bam_sort, bl, bi);
		continue;
	    }

	    if (bi->keep_qual)
		memcpy(bam_get_qual(bi->b), bam_get_qual(plp[i].b), plp[i].b->core.l_qseq);

	    // Sliding window to not over-egg the pudding.
	    // If a region of a read is obviously duff, then even if it
	    // matches we should be cautious of over boosting the qual.
	    // Instead we use binning instead, allowing lowering.
	    if (0) {
		int qa = 0, x;
		bam1_t *b = plp[i].b;
		uint8_t *qual = bam_get_qual(b2);
		uint8_t *qual_orig = bam_get_qual(b);
		int WL = 10;
		for (x = 0; x < WL && x < b->core.l_qseq; x++) {
		    qa += qual_orig[x] < 15;
		}
		while (x < b->core.l_qseq-1) {
		    if (qa >= 6) {
			//memcpy(&qual[x-WL], &qual_orig[x-WL], WL);
			int j;
			for (j = x-WL; j < x; j++)
			    qual[j] = bin2[qual_orig[j]];
		    }
		    qa -= qual_orig[x-WL] < 15;
		    qa += qual_orig[x++] < 15;
		}
	    }

	    uint64_t id = bi->id;
	    bam_sorted_item *next = RB_NEXT(bam_sort, bl, bi);
	    remove_bam_list(bl, bi);
	    bi = next;

	    // Note: this may reorder seqs that start at the same coord,
	    // so we give it the read-id to preserve the order.
	    insert_bam_list_id(b_hist, b2, id);
	}

	// Flush history (preserving sort order).
	if (flush_bam_list(&cd, p, b_hist, INT_MAX, left_most, out, header) < 0)
	    return -1;
    }

    // Handle any in-flight reads that haven't yet finished as pileup
    // was called with a range and we've terminated the pileup iterator.
    bam_sorted_item *bi = RB_MIN(bam_sort, bl);
    while (bi) {
	bam_sorted_item *next = RB_NEXT(bam_sort, bl, bi);
	insert_bam_list_id(b_hist, bi->b, bi->id);
	remove_bam_list(bl, bi);
	bi = next;
    }

    if (flush_bam_list(&cd, p, b_hist, INT_MAX, INT_MAX, out, header) < 0)
	return -1;

    // Handle trailing unmapped reads
    if (cd.b_unmap) {
	int next = 0;
	do {
	    purge_tags(p, cd.b_unmap);
	    cd.count_out++;
	    if (p->pblock)
		pblock(cd.b_unmap, p->pblock);
	    if (sam_write1(out, header, cd.b_unmap) < 0)
		return -1;
	    next = (cd.iter
		    ? sam_itr_next(cd.fp, cd.iter, cd.b_unmap)
		    : sam_read1(cd.fp, cd.header, cd.b_unmap)) >= 0;
	    if (next) cd.count_in++;
	} while (next);

	bam_destroy1(cd.b_unmap);
    }

    bam_plp_destroy(p_iter);
    bam_sorted_list_destroy(bl);
    bam_sorted_list_destroy(b_hist);

    if (cd.count_in != cd.count_out) {
	fprintf(stderr, "ERROR: lost a read?\n");
	fprintf(stderr, "Read  %"PRId64" reads\n",   cd.count_in);
	fprintf(stderr, "Wrote %"PRId64" reads\n\n", cd.count_out);
	return 1;
    }

    return 0;
}

int parse_aux_list(auxhash_t *h, char *optarg) {
    if (!*h)
        *h = kh_init(aux_exists);

    while (strlen(optarg) >= 2) {
        int x = optarg[0]<<8 | optarg[1];
        int ret = 0;
        kh_put(aux_exists, *h, x, &ret);

        optarg += 2;
        if (*optarg == ',') // allow white-space too for easy `cat file`?
            optarg++;
        else if (*optarg != 0)
            break;
    }

    if (strlen(optarg) != 0) {
        fprintf(stderr, "main_samview: Error parsing option, "
                "auxiliary tags should be exactly two characters long.\n");
        return -1;
    }

    return 0;
}

void usage(FILE *fp) {
    fprintf(fp, "Crumble version %s\n\n", CRUMBLE_VERSION);
    fprintf(fp, "Usage: crumble [options] in-file out-file\n");
    fprintf(fp, "\nOptions:\n"
	    "-I fmt(,opt...)   Input format and format-options [auto].\n"
	    "-O fmt(,opt...)   Output format and format-options [SAM].\n");
    fprintf(fp, "-v                Increase verbosity\n");
    fprintf(fp, "-z                Do not add an @PG SAM header line\n");
    fprintf(fp,
"-c qual_cutoff    In highly confident regions, quality values above/below\n"
"-l qual_lower         'qual_cutoff' [%d] are quantised to 'qual_lower' [%d]\n"
"-u qual_upper         and 'qual_upper' [%d] based on agreement to consensus.\n",
	    QM, QL, QH);
    fprintf(fp, "-S                Quantise qualities (with -[clu] options) in soft-clips too.\n");
    fprintf(fp, "-m min_mqual      Keep qualities for seqs with mapping quality <= mqual [%d].\n", MIN_MQUAL);
    fprintf(fp, "-L bool           Whether mismatching bases can have qualities lowered [%d]\n", REDUCE_QUAL);
    fprintf(fp, "-B                If set, replace quals in good regions with low/high [%sset]\n", BINARY_QUAL ? "" : "un");
    fprintf(fp, "-i STR_mul,add    Adjust indel size by (STR_size+add)*mul [%.1f,%d]\n", I_STR_MUL, I_STR_ADD);
    fprintf(fp, "-s STR_mul,add    Adjust SNP size by (STR_size+add)*mul [%.1f,%d]\n", S_STR_MUL, S_STR_ADD);
    fprintf(fp, "-r region         Limit input to region chr:pos(-pos) []\n");
    fprintf(fp, "-t tag_list       Comma separated list of aux tags to keep []\n");
    fprintf(fp, "-T tag_list       Comma separated list of aux tags to discard []\n");
    fprintf(fp, "-b out.bed        Output suspicious regions to out.bed []\n");
    fprintf(fp, "-P float          Keep qual if depth locally >= [%.1f] times deeper than expected\n", OVER_DEPTH);
    fprintf(fp, "\n(Preserving whole read qualities; values are fractions of read coverage)\n");
    fprintf(fp, "-C float          Keep if >= [%.2f] reads have soft-clipping\n", CLIP_PERC);
    fprintf(fp, "-M float          Keep if >= [%.2f] reads have low mapping quality\n", LOW_MQUAL_PERC);
    fprintf(fp, "-Z float          Keep if >= [%.2f] indel sizes do not fit bi-modal dist.\n", INS_LEN_PERC);
    fprintf(fp, "-V float          Keep if <  [%.2f] reads span indel\n", INDEL_OVERLAP_PERC);
    fprintf(fp, "\n(Calling while ignoring mapping quality)\n");
    fprintf(fp, "-q int            Minimum snp call confidence [%d]\n", MIN_QUAL_A);
    fprintf(fp, "-d int            Minimum indel call confidence [%d]\n", MIN_INDEL_A);
    fprintf(fp, "-x float          Minimum discrepancy score [%.1f]\n", MIN_DISCREP_A);
    fprintf(fp, "\n(Calling with use of mapping quality)\n");
    fprintf(fp, "-Q int            Minimum snp call confidence [%d]\n", MIN_QUAL_B);
    fprintf(fp, "-D int            Minimum indel call confidence [%d]\n", MIN_INDEL_B);
    fprintf(fp, "-X float          Minimum discrepancy score [%.1f]\n", MIN_DISCREP_B);
    fprintf(fp, "\n(Horizontal quality smoothing via P-block)\n");
    fprintf(fp, "-p int            P-block algorithm; quality values +/- 'int' [0]\n");
    fprintf(fp, "\n(BD and BI aux tag binary-binning; off by default)\n");
    fprintf(fp, "-f qual_cutoff    Quantise BD:Z: tags to two values (or one if both equal).\n");
    fprintf(fp, "-g qual_upper       If >= 'qual_cutoff' [0] replace by 'qual_upper' [0]\n");
    fprintf(fp, "-e qual_lower       otherwise replace by 'qual_lower' [0].\n");
    fprintf(fp, "-F qual_cutoff    Quantise BI:Z: tags to two values (or one if both equal).\n");
    fprintf(fp, "-G qual_upper       If >= 'qual_cutoff' [0] replace by 'qual_upper' [0]\n");
    fprintf(fp, "-E qual_lower       otherwise replace by 'qual_lower' [0].\n");
    fprintf(fp, "\n(Standard compression levels combining the above.)\n");
    fprintf(fp, "-1                Synonym for -s1.0,5 -i2.0,1 -m5\n");
    fprintf(fp, "-3                Synonym for -s1.0,0\n");
    fprintf(fp, "-5                Synonym for (defaults)\n");
    fprintf(fp, "-7                Synonym for -P 999 -C 1 -M 1 -Z 1 -V 0\n");
    fprintf(fp, "-9                Synonym for -Q70 -D125 -X1.5 -P 999 -C 1 -M 1 -Z 1 -V 0\n");
    fprintf(fp, "\n");
    fprintf(fp,
"Standard htslib format options apply.  So to create a CRAM file with lossy\n\
template names enabled and a larger number of sequences per slice, try:\n\
\n\
    crumble -O cram,lossy_names,seqs_per_slice=100000\n\
\n\
The lossy quality encoding works by running two distinct heterozygous consensus\n\
calling algorithms; with and without the use of mapping qualities.  Use -q 0\n\
or -Q 0 to disable one of these if only the other is needed.  When operating,\n\
any sufficiently high quality SNP (above -q / -Q) with have the qualities for\n\
the bases adjusted to 'qual_lower' or 'qual_upper'.  Similarly for any high\n\
quality indel.  An lower quality indel will causes neighbouring bases for\n\
all sequences at that site to be kept, for the region as large as the indel\n\
plus an extension along any short tandem repeats (STR), multiplied by \n\
'indel_mult' plus an additional 'STR_add'.\n");
}

int main(int argc, char **argv) {
    samFile *in, *out = NULL;
    htsFormat in_fmt = {0};
    htsFormat out_fmt = {0};
    bam_hdr_t *header;
    hts_itr_t *h_iter = NULL;
    int opt;

    cram_lossy_params params = {
	.reduce_qual   = REDUCE_QUAL,       // -L
	.binary_qual   = BINARY_QUAL,       // -B
	.iSTR_mul      = I_STR_MUL,         // -i
	.iSTR_add      = I_STR_ADD,         // -i
	.sSTR_mul      = S_STR_MUL,         // -s
	.sSTR_add      = S_STR_ADD,         // -s
	.qlow          = QL,                // -l
	.qcutoff       = QM,		    // -c
	.qhigh         = QH,		    // -u
	.min_mqual     = MIN_MQUAL,	    // -m
	.min_qual_A    = MIN_QUAL_A,	    // -q
	.min_indel_A   = MIN_INDEL_A,	    // -d
	.min_discrep_A = MIN_DISCREP_A,	    // -x
	.min_qual_B    = MIN_QUAL_B,	    // -Q
	.min_indel_B   = MIN_INDEL_B,	    // -D
	.min_discrep_B = MIN_DISCREP_B,     // -X
	.aux_whitelist = NULL,              // -t
	.aux_blacklist = NULL,              // -T
	.region        = NULL,              // -r
	.bed_fp        = NULL,              // -b
	.clip_perc     = CLIP_PERC,         // -C
	.low_mqual_perc= LOW_MQUAL_PERC,    // -M
	.ins_len_perc  = INS_LEN_PERC,      // -Z
	.over_depth    = OVER_DEPTH,        // -P
	.indel_ov_perc = INDEL_OVERLAP_PERC,// -V
	.verbose       = 0,                 // -v
	.pblock        = 0,                 // -p
	.BD_low        = 0,                 // -e
	.BD_mid        = 0,                 // -f
	.BD_high       = 0,                 // -g
	.BI_low        = 0,                 // -E
	.BI_mid        = 0,                 // -F
	.BI_high       = 0,                 // -G
	.softclip      = 0,                 // -S
	.noPG          = 0,                 // -z
    };

    //  ........  ..  ....... . .
    // abcdefghijklmnopqrstuvwxyz
    //  ...... .  .. ... .. . . .
    // ABCDEFGHIJKLMNOPQRSTUVWXYZ

    while ((opt = getopt(argc, argv, "I:O:q:d:x:Q:D:X:m:l:u:c:i:L:Bs:t:T:hr:b:vC:M:Z:P:V:p:e:f:g:E:F:G:S13579z")) != -1) {
	switch (opt) {
	case 'I':
	    hts_parse_format(&in_fmt, optarg);
	    break;

	case 'O':
	    hts_parse_format(&out_fmt, optarg);
	    break;

	case 'q':
	    params.min_qual_A = atoi(optarg);
	    break;
	case 'd':
	    params.min_indel_A = atoi(optarg);
	    break;
	case 'x':
	    params.min_discrep_A = atof(optarg);
	    break;

	case 'Q':
	    params.min_qual_B = atoi(optarg);
	    break;
	case 'D':
	    params.min_indel_B = atoi(optarg);
	    break;
	case 'X':
	    params.min_discrep_B = atof(optarg);
	    break;

	case 'm':
	    params.min_mqual = atoi(optarg);
	    break;

	case 'l':
	    params.qlow = atoi(optarg);
	    break;
	case 'u':
	    params.qhigh = atoi(optarg);
	    break;
	case 'c':
	    params.qcutoff = atoi(optarg);
	    break;

	case 'i':
	    params.iSTR_mul = atof(optarg);
	    if (strchr(optarg,','))
		params.iSTR_add = atoi(strchr(optarg,',')+1);
	    break;

	case 's':
	    params.sSTR_mul = atof(optarg);
	    if (strchr(optarg,','))
		params.sSTR_add = atoi(strchr(optarg,',')+1);
	    break;

	case 'L':
	    params.reduce_qual = atoi(optarg);
	    break;

	case 'B':
	    params.binary_qual = 1;
	    break;

	case 'r':
	    params.region = optarg;
	    break;

	case 't':
            if (parse_aux_list(&params.aux_whitelist, optarg)) {
		usage(stderr);
                return 1;
	    }
            break;

	case 'T':
            if (parse_aux_list(&params.aux_blacklist, optarg)) {
		usage(stderr);
                return 1;
	    }
            break;

	case 'b':
	    if ((params.bed_fp = fopen(optarg, "w")) == NULL) {
		perror(optarg);
		return 1;
	    }
	    break;

	case 'C':
	    params.clip_perc = atof(optarg);
	    break;

	case 'M':
	    params.low_mqual_perc = atof(optarg);
	    break;

	case 'Z':
	    params.ins_len_perc = atof(optarg);
	    break;

	case 'P':
	    params.over_depth = atof(optarg);
	    break;

	case 'V':
	    params.indel_ov_perc = atof(optarg);
	    break;

	case 'p':
	    params.pblock = atoi(optarg);
	    break;

	case 'e':
	    params.BD_low = atoi(optarg)+33;
	    break;

	case 'f':
	    params.BD_mid = atoi(optarg)+33;
	    break;

	case 'g':
	    params.BD_high = atoi(optarg)+33;
	    break;

	case 'E':
	    params.BI_low = atoi(optarg)+33;
	    break;

	case 'F':
	    params.BI_mid = atoi(optarg)+33;
	    break;

	case 'G':
	    params.BI_high = atoi(optarg)+33;
	    break;

	case '9':
	    // Most aggressive compression
	    params.min_qual_B = 70;
	    params.min_indel_B = 125;
	    params.min_discrep_B = 1.5;
	    params.low_mqual_perc = 1.0;
	    params.ins_len_perc = 1.0;
	    params.indel_ov_perc = 0;
	    params.over_depth = 999;
	    params.iSTR_mul = 1.0;
	    break;

	case '7':
	    params.low_mqual_perc = 1.0;
	    params.ins_len_perc = 1.0;
	    params.indel_ov_perc = 0;
	    params.over_depth = 999;
	    break;

	case '5':
	    break;

	case '3':
	    params.sSTR_mul = 1.0;
	    params.sSTR_add = 0;
	    break;

	case '1':
	    // Most conservative compression
	    params.sSTR_mul = 1.0;
	    params.sSTR_add = 5;
	    params.iSTR_mul = 2.0;
	    params.iSTR_add = 1;
	    params.min_mqual = 5;
	    break;

	case 'S':
	    params.softclip = 1;
	    break;

	case 'z':
	    params.noPG = 1;
	    break;

	case 'v':
	    params.verbose++;
	    break;

	case 'h':
	    usage(stdout);
	    return 1;

	default: /* ? */
	    usage(stderr);
	    return 1;
	}
    }

    if (params.verbose) {
	printf("--- Crumble v%s: parameters ---\n", CRUMBLE_VERSION);
	printf("reduce qual:   %s\n",     params.reduce_qual ? "yes" : "no");
	printf("indel STR mul: %.2f\n",   params.iSTR_mul);
	printf("indel STR add: %d\n",     params.iSTR_add);
	printf("SNP   STR mul: %.2f\n",   params.sSTR_mul);
	printf("SNP   STR add: %d\n",     params.sSTR_add);
	if (params.binary_qual) {
	    printf("Qual low  1..%d -> %d\n", params.qcutoff-1, params.qlow);
	    printf("Qual high %d..  -> %d\n", params.qcutoff, params.qhigh);
	} else {
	    printf("Qual low  %d, used for discrepant bases in high conf call\n", params.qlow);
	    printf("Qual high %d, used for matching bases in high conf call\n", params.qhigh);
	}
	printf("Keep if mqual <= %d\n",   params.min_mqual);
	if (params.min_qual_A) {
	    printf("Calls without mqual, keep qual if:\n");
	    printf("  SNP < %d,  indel < %d,  discrep > %.2f\n",
		   params.min_qual_A, params.min_indel_A, params.min_discrep_A);
	} else {
	    printf("Calls without mqual: disabled.\n");
	}
	if (params.min_qual_B) {
	    printf("Calls with mqual, keep qual if:\n");
	    printf("  SNP < %d,  indel < %d,  discrep > %.2f\n",
		   params.min_qual_B, params.min_indel_B, params.min_discrep_B);
	} else {
	    printf("Calls with mqual: disabled.\n");
	}
    }

    init_bins(&params);

    char *fnin = NULL;
    if ( optind>=argc ) {
        if ( !isatty(STDIN_FILENO) ) fnin = "-";  // reading from stdin
        else { usage(stdout); return 1; }
    }
    else fnin = argv[optind++];

    if (!(in = sam_open_format(fnin, "r", &in_fmt))) {
	perror(argv[optind]);
	return 1;
    }

    char mode[5] = "w";
    char *fnout = optind < argc ? argv[optind++] : "-";
    sam_open_mode(mode+1, fnout, NULL);

    if (!(out = sam_open_format(fnout, mode, &out_fmt))) {
	perror("(stdout)");
	return 1;
    }

    if (!(header = sam_hdr_read(in))) {
	fprintf(stderr, "Failed to read file header\n");
	return 1;
    }

    if (!params.noPG) {
	// Abuse the internal CRAM htslib header parser; to become public in due course.
	SAM_hdr *sh = sam_hdr_parse_(header->text, header->l_text);
	if (!sh)
	    return 1;

	char *arg_list = stringify_argv(argc, argv);

	if (sam_hdr_add_PG(sh, "crumble",
			   "VN", CRUMBLE_VERSION,
			   arg_list ? "CL" : NULL,
			   arg_list ? arg_list : NULL,
			   NULL) != 0)
	    return 1;

	free(header->text);
	if (!(header->text = strdup(sam_hdr_str(sh))))
	    return 1;
	header->l_text = sam_hdr_length(sh);
	sam_hdr_free(sh);
	free(arg_list);
    }

    if (out && sam_hdr_write(out, header) != 0) {
	fprintf(stderr, "Failed to write file header\n");
	return 1;
    }

    if (params.region) {
        hts_idx_t *idx = sam_index_load(in, fnin);
    	h_iter = idx ? sam_itr_querys(idx, header, params.region) : NULL;
    	if (!h_iter || !idx) {
    	    fprintf(stderr, "Failed to load index and/or parse iterator.\n");
    	    return 1;
    	}
    	hts_idx_destroy(idx);
    }

    if (transcode(&params, in, out, header, h_iter) != 0) {
	fprintf(stderr, "Error while reducing file\n");
	return 1;
    }

    bam_hdr_destroy(header);
    if (h_iter) hts_itr_destroy(h_iter);

    if (sam_close(in) != 0) {
	fprintf(stderr, "Error while closing input fd\n");
	return 1;
    }

    if (out && sam_close(out) != 0) {
	fprintf(stderr, "Error while closing output fd\n");
	return 1;
    }

    if (params.aux_whitelist)
	kh_destroy(aux_exists, params.aux_whitelist);

    if (params.aux_blacklist)
	kh_destroy(aux_exists, params.aux_blacklist);

    if (params.verbose) {
	fprintf(stderr, "A/B Diff         = %d\n", count_diff);
	fprintf(stderr, "A/B Indel        = %d / %d\n", count_indel_qual, count_indel); 
	fprintf(stderr, "A:  Het          = %d / %d\n", count_het_qual_A, count_het_A);
	fprintf(stderr, "A:  Hom          = %d / %d\n", count_hom_qual_A, count_hom_A);
	fprintf(stderr, "A:  Discrep      = %d\n", count_discrep_A);
	fprintf(stderr, "B:  Het          = %d / %d\n", count_het_qual_B, count_het_B);
	fprintf(stderr, "B:  Hom          = %d / %d\n", count_hom_qual_B, count_hom_B);
	fprintf(stderr, "B:  Discrep      = %d\n\n", count_discrep_B);
	fprintf(stderr, "Columns          = %"PRId64"\n", count_columns);
	fprintf(stderr, "Low_mqual_perc   = %"PRId64"\n", count_low_mqual_perc);
	fprintf(stderr, "Clip_perc        = %"PRId64"\n", count_clip_perc);
	fprintf(stderr, "Ins_len_perc     = %"PRId64"\n", count_ins_len_perc);
	fprintf(stderr, "indel_ov_perc    = %"PRId64"\n", count_indel_ov_perc);
	fprintf(stderr, "count_over_depth = %"PRId64"\n", count_over_depth);
    }

    if (params.bed_fp)
	fclose(params.bed_fp);

    return 0;
}
