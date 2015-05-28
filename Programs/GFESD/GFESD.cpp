// Updated on 05/25/15 

/* 

Program GFESD.cpp to:

	1) generate site-specific read data for random individuals, with given genotype frequencies and error rates.
	2) use ML to estimate the genotype frequencies and error rate, in retrospect.

Data are for multiple diploid individuals from a population where Hardy-Weinberg equilibrium may be violated.
Coverage per site is assumed to be Poisson distributed among individuals.

In the simulation for generating read data, the major allele is designated 1, and the minor allele 2, with 3 and 4 being errors to the other 
	two bases.

*/


/* ******************************************************************************************* */

#include	<stdio.h>
#include 	<math.h>
#include 	<sys/types.h>
#include	<stdlib.h>
#include	<time.h>

/* ******************************************************************************************* */




/* ALL OF THE FOLLOWING MATERIAL, UP TO THE MAIN PROGRAM, IS USED FOR GENERATING RANDOM NUMBERS */

extern long ignbin(long n,double pp);


/* Definitions for the binomial generator. */

#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))


/* Definitions for the random number generator. */

#define TRUE -1
#define FALSE 0
#define boolean int


static double u[98], c, cd, cm;
static int i97, j97;
static boolean test = FALSE;


double ranmar(void)
/*
C This is the random number generator proposed by George Marsaglia in 
C Florida State University Report: FSU-SCRI-87-50
C It was slightly modified by F. James to produce an array of pseudorandom
C numbers.
*/
{
        double uni;
        
        if (test==FALSE) {
                puts("Call the init routine rmarin() before calling ranmar().");
                exit(2);
        }
	uni = u[i97] - u[j97];
	if (uni < 0.0) uni += 1.0;
	u[i97] = uni;
	i97--;
	if (i97==0) i97 = 97;
	j97--;
	if (j97==0) j97 = 97;
	c -= cd;
	if (c<0.0) c += cm;
	uni -= c;
	if (uni<0.0) uni += 1.0;
	return uni;
}




/* Seed for the random number generator. */

void rmarin(int ij,int kl) 
{
/*
C This is the initialization routine for the random number generator RANMAR()
C NOTE: The seed variables can have values between:    0 <= IJ <= 31328
C                                                      0 <= KL <= 30081
C The random number sequences created by these two seeds are of sufficient 
C length to complete an entire calculation with. For example, if sveral 
C different groups are working on different parts of the same calculation,
C each group could be assigned its own IJ seed. This would leave each group
C with 30000 choices for the second seed. That is to say, this random 
C number generator can create 900 million different subsequences -- with 
C each subsequence having a length of approximately 10^30.
C 
C Use IJ = 1802 & KL = 9373 to test the random number generator. The
C subroutine RANMAR should be used to generate 20000 random numbers.
C Then display the next six random numbers generated multiplied by 4096*4096
C If the random number generator is working properly, the random numbers
C should be:
C           6533892.0  14220222.0  7275067.0
C           6172232.0  8354498.0   10633180.0
*/
        int i, j, k, l, ii, jj, m;
        double s, t;
        
        if (ij<0 || ij>31328 || kl<0 || kl>30081) {
                puts("The first random number seed must have a value between 0 and 31328.");
                puts("The second seed must have a value between 0 and 30081.");
                exit(1);
        }
        
        i = (ij/177)%177 + 2;
        j = ij%177 + 2;
        k = (kl/169)%178 + 1;
        l = kl%169;
        
        for (ii=1; ii<=97; ii++) {
                s = 0.0;
                t = 0.5;
                for (jj=1; jj<=24; jj++) {
                        m = (((i*j)%179)*k) % 179;
                        i = j;
                        j = k;
                        k = m;
                        l = (53*l + 1) % 169;
                        if ((l*m)%64 >= 32) s += t;
                        t *= 0.5;
                }
                u[ii] = s;
        }
        
        c = 362436.0 / 16777216.0;
        cd = 7654321.0 / 16777216.0;
        cm = 16777213.0 / 16777216.0;
        
        i97 = 97;
        j97 = 33;
        
        test = TRUE;
}




/* ****************************************************************************************** */

/* binomial random number generator */

long ignbin(long n,double pp)
/*
**********************************************************************
     long ignbin(long n,double pp)
                    GENerate BINomial random deviate
                              Function
     Generates a single random deviate from a binomial
     distribution whose number of trials is N and whose
     probability of an event in each trial is P.
                              Arguments
     n  --> The number of trials in the binomial distribution
            from which a random deviate is to be generated.
     p  --> The probability of an event in each trial of the
            binomial distribution from which a random deviate
            is to be generated.
     ignbin <-- A random deviate yielding the number of events
                from N independent trials, each of which has
                a probability of event P.
                              Method
     This is algorithm BTPE from:
         Kachitvichyanukul, V. and Schmeiser, B. W.
         Binomial Random Variate Generation.
         Communications of the ACM, 31, 2
         (February, 1988) 216.
**********************************************************************
     SUBROUTINE BTPEC(N,PP,ISEED,JX)
     BINOMIAL RANDOM VARIATE GENERATOR
     MEAN .LT. 30 -- INVERSE CDF
       MEAN .GE. 30 -- ALGORITHM BTPE:  ACCEPTANCE-REJECTION VIA
       FOUR REGION COMPOSITION.  THE FOUR REGIONS ARE A TRIANGLE
       (SYMMETRIC IN THE CENTER), A PAIR OF PARALLELOGRAMS (ABOVE
       THE TRIANGLE), AND EXPONENTIAL LEFT AND RIGHT TAILS.
     BTPE REFERS TO BINOMIAL-TRIANGLE-PARALLELOGRAM-EXPONENTIAL.
     BTPEC REFERS TO BTPE AND "COMBINED."  THUS BTPE IS THE
       RESEARCH AND BTPEC IS THE IMPLEMENTATION OF A COMPLETE
       USABLE ALGORITHM.
     REFERENCE:  VORATAS KACHITVICHYANUKUL AND BRUCE SCHMEISER,
       "BINOMIAL RANDOM VARIATE GENERATION,"
       COMMUNICATIONS OF THE ACM, FORTHCOMING
     WRITTEN:  SEPTEMBER 1980.
       LAST REVISED:  MAY 1985, JULY 1987
     REQUIRED SUBPROGRAM:  RAND() -- A UNIFORM (0,1) RANDOM NUMBER
                           GENERATOR
     ARGUMENTS
       N : NUMBER OF BERNOULLI TRIALS            (INPUT)
       PP : PROBABILITY OF SUCCESS IN EACH TRIAL (INPUT)
       ISEED:  RANDOM NUMBER SEED                (INPUT AND OUTPUT)
       JX:  RANDOMLY GENERATED OBSERVATION       (OUTPUT)
     VARIABLES
       PSAVE: VALUE OF PP FROM THE LAST CALL TO BTPEC
       NSAVE: VALUE OF N FROM THE LAST CALL TO BTPEC
       XNP:  VALUE OF THE MEAN FROM THE LAST CALL TO BTPEC
       P: PROBABILITY USED IN THE GENERATION PHASE OF BTPEC
       FFM: TEMPORARY VARIABLE EQUAL TO XNP + P
       M:  INTEGER VALUE OF THE CURRENT MODE
       FM:  FLOATING POINT VALUE OF THE CURRENT MODE
       XNPQ: TEMPORARY VARIABLE USED IN SETUP AND SQUEEZING STEPS
       P1:  AREA OF THE TRIANGLE
       C:  HEIGHT OF THE PARALLELOGRAMS
       XM:  CENTER OF THE TRIANGLE
       XL:  LEFT END OF THE TRIANGLE
       XR:  RIGHT END OF THE TRIANGLE
       AL:  TEMPORARY VARIABLE
       XLL:  RATE FOR THE LEFT EXPONENTIAL TAIL
       XLR:  RATE FOR THE RIGHT EXPONENTIAL TAIL
       P2:  AREA OF THE PARALLELOGRAMS
       P3:  AREA OF THE LEFT EXPONENTIAL TAIL
       P4:  AREA OF THE RIGHT EXPONENTIAL TAIL
       U:  A U(0,P4) RANDOM VARIATE USED FIRST TO SELECT ONE OF THE
           FOUR REGIONS AND THEN CONDITIONALLY TO GENERATE A VALUE
           FROM THE REGION
       V:  A U(0,1) RANDOM NUMBER USED TO GENERATE THE RANDOM VALUE
           (REGION 1) OR TRANSFORMED INTO THE VARIATE TO ACCEPT OR
           REJECT THE CANDIDATE VALUE
       IX:  INTEGER CANDIDATE VALUE
       X:  PRELIMINARY CONTINUOUS CANDIDATE VALUE IN REGION 2 LOGIC
           AND A FLOATING POINT IX IN THE ACCEPT/REJECT LOGIC
       K:  ABSOLUTE VALUE OF (IX-M)
       F:  THE HEIGHT OF THE SCALED DENSITY FUNCTION USED IN THE
           ACCEPT/REJECT DECISION WHEN BOTH M AND IX ARE SMALL
           ALSO USED IN THE INVERSE TRANSFORMATION
       R: THE RATIO P/Q
       G: CONSTANT USED IN CALCULATION OF PROBABILITY
       MP:  MODE PLUS ONE, THE LOWER INDEX FOR EXPLICIT CALCULATION
            OF F WHEN IX IS GREATER THAN M
       IX1:  CANDIDATE VALUE PLUS ONE, THE LOWER INDEX FOR EXPLICIT
             CALCULATION OF F WHEN IX IS LESS THAN M
       I:  INDEX FOR EXPLICIT CALCULATION OF F FOR BTPE
       AMAXP: MAXIMUM ERROR OF THE LOGARITHM OF NORMAL BOUND
       YNORM: LOGARITHM OF NORMAL BOUND
       ALV:  NATURAL LOGARITHM OF THE ACCEPT/REJECT VARIATE V
       X1,F1,Z,W,Z2,X2,F2, AND W2 ARE TEMPORARY VARIABLES TO BE
       USED IN THE FINAL ACCEPT/REJECT TEST
       QN: PROBABILITY OF NO SUCCESS IN N TRIALS
     REMARK
       IX AND JX COULD LOGICALLY BE THE SAME VARIABLE, WHICH WOULD
       SAVE A MEMORY POSITION AND A LINE OF CODE.  HOWEVER, SOME
       COMPILERS (E.G.,CDC MNF) OPTIMIZE BETTER WHEN THE ARGUMENTS
       ARE NOT INVOLVED.
     ISEED NEEDS TO BE DOUBLE PRECISION IF THE IMSL ROUTINE
     GGUBFS IS USED TO GENERATE UNIFORM RANDOM NUMBER, OTHERWISE
     TYPE OF ISEED SHOULD BE DICTATED BY THE UNIFORM GENERATOR
**********************************************************************
*****DETERMINE APPROPRIATE ALGORITHM AND WHETHER SETUP IS NECESSARY
*/

{
static float psave = -1.0;
static long nsave = -1;
static long ignbin,i,ix,ix1,k,m,mp,T1;
static double al,alv,amaxp,c,f,f1,f2,ffm,fm,g,p,p1,p2,p3,p4,q,qn,r,u,v,w,w2,x,x1,
    x2,xl,xll,xlr,xm,xnp,xnpq,xr,ynorm,z,z2;

    if(pp != psave) goto S10;
    if(n != nsave) goto S20;
    if(xnp < 30.0) goto S150;
    goto S30;
S10:
/*
*****SETUP, PERFORM ONLY WHEN PARAMETERS CHANGE
*/
    psave = pp;
    p = min(psave,1.0-psave);
    q = 1.0-p;
S20:
    xnp = n*p;
    nsave = n;
    if(xnp < 30.0) goto S140;
    ffm = xnp+p;
    m = (long)ffm;
    fm = m;
    xnpq = xnp*q;
    p1 = (long) (2.195*sqrt(xnpq)-4.6*q)+0.5;
    xm = fm+0.5;
    xl = xm-p1;
    xr = xm+p1;
    c = 0.134+20.5/(15.3+fm);
    al = (ffm-xl)/(ffm-xl*p);
    xll = al*(1.0+0.5*al);
    al = (xr-ffm)/(xr*q);
    xlr = al*(1.0+0.5*al);
    p2 = p1*(1.0+c+c);
    p3 = p2+c/xll;
    p4 = p3+c/xlr;
S30:
/*
*****GENERATE VARIATE
*/
    u = ranmar()*p4;
    v = ranmar();
/*
     TRIANGULAR REGION
*/
    if(u > p1) goto S40;
    ix = (long)(xm-p1*v+u);
    goto S170;
S40:
/*
     PARALLELOGRAM REGION
*/
    if(u > p2) goto S50;
    x = xl+(u-p1)/c;
    v = v*c+1.0-ABS(xm-x)/p1;
    if(v > 1.0 || v <= 0.0) goto S30;
    ix = (long)x;
    goto S70;
S50:
/*
     LEFT TAIL
*/
    if(u > p3) goto S60;
    ix = (long)(xl+log(v)/xll);
    if(ix < 0) goto S30;
    v *= ((u-p2)*xll);
    goto S70;
S60:
/*
     RIGHT TAIL
*/
    ix = (long)(xr-log(v)/xlr);
    if(ix > n) goto S30;
    v *= ((u-p3)*xlr);
S70:
/*
*****DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST
*/
    k = ABS(ix-m);
    if(k > 20 && k < xnpq/2-1) goto S130;
/*
     EXPLICIT EVALUATION
*/
    f = 1.0;
    r = p/q;
    g = (n+1)*r;
    T1 = m-ix;
    if(T1 < 0) goto S80;
    else if(T1 == 0) goto S120;
    else  goto S100;
S80:
    mp = m+1;
    for(i=mp; i<=ix; i++) f *= (g/i-r);
    goto S120;
S100:
    ix1 = ix+1;
    for(i=ix1; i<=m; i++) f /= (g/i-r);
S120:
    if(v <= f) goto S170;
    goto S30;
S130:
/*
     SQUEEZING USING UPPER AND LOWER BOUNDS ON ALOG(F(X))
*/
    amaxp = k/xnpq*((k*(k/3.0+0.625)+0.1666666666666)/xnpq+0.5);
    ynorm = -(k*k/(2.0*xnpq));
    alv = log(v);
    if(alv < ynorm-amaxp) goto S170;
    if(alv > ynorm+amaxp) goto S30;
/*
     STIRLING'S FORMULA TO MACHINE ACCURACY FOR
     THE FINAL ACCEPTANCE/REJECTION TEST
*/
    x1 = ix+1.0;
    f1 = fm+1.0;
    z = n+1.0-fm;
    w = n-ix+1.0;
    z2 = z*z;
    x2 = x1*x1;
    f2 = f1*f1;
    w2 = w*w;
    if(alv <= xm*log(f1/x1)+(n-m+0.5)*log(z/w)+(ix-m)*log(w*p/(x1*q))+(13860.0-
      (462.0-(132.0-(99.0-140.0/f2)/f2)/f2)/f2)/f1/166320.0+(13860.0-(462.0-
      (132.0-(99.0-140.0/z2)/z2)/z2)/z2)/z/166320.0+(13860.0-(462.0-(132.0-
      (99.0-140.0/x2)/x2)/x2)/x2)/x1/166320.0+(13860.0-(462.0-(132.0-(99.0
      -140.0/w2)/w2)/w2)/w2)/w/166320.0) goto S170;
    goto S30;
S140:
/*
     INVERSE CDF LOGIC FOR MEAN LESS THAN 30
*/
    qn = pow(q,(double)n);
    r = p/q;
    g = r*(n+1);
S150:
    ix = 0;
    f = qn;
    u = ranmar();
S160:
    if(u < f) goto S170;
    if(ix > 110) goto S150;
    u -= f;
    ix += 1;
    f *= (g/ix-r);
    goto S160;
S170:
    if(psave > 0.5) ix = n-ix;
    ignbin = ix;
    return ignbin;
}


/* ****************************************************************************************** */





// point to the output file 

FILE *stream;

int main(int argc, char *argv[])

{

int nsample, coverage, niters;
float epsilon;
double P, Q, D_A, f;
char out_file_name[100];

// Read the specified setting
if (argc != 8) {	// A check to make sure the user has passed the right number of arguments.
	fprintf(stderr, "Usage : %s N e n reps P Q out_file_name\n", argv[0]); 
	exit(1);
}
sscanf(argv[1], "%d", &nsample);
sscanf(argv[2], "%f", &epsilon);
sscanf(argv[3], "%d", &coverage);
sscanf(argv[4], "%d", &niters);
sscanf(argv[5], "%lf", &P);
sscanf(argv[6], "%lf", &Q);
sscanf(argv[7], "%s", out_file_name);

// Calculate the true values of the disequilibrium coefficient D_A and inbreeding coefficient f
D_A = ( (1.0+P-Q)*(1.0-P+Q) )/4.0 - (1.0 - P - Q)/2.0;
if (P < 1.0 - 1.0e-8) {
	f = 1.0 - (1.0 - P - Q)/( 0.5*(1.0+P-Q)*(1.0-P+Q) );
} else {
	f = 0.0;
}

/* NOTE: The seed variables can have values between:    0 <= IJ <= 31328 */
/*                                                      0 <= KL <= 30081 */
/* Default random seeds for the ranmar() random number generator:        */

int ij = 1802;
int kl = 9373;

/* Uncomment the next line for system clock seeding of ranmar(): */
 ij = time ((long *) 0) % 31329; 

/* Initialization of ranmar() */
/* i.e., must always call rmarin(ij,kl) before first usage of ranmar */

rmarin(ij,kl);







/* ************************************************************************************************************ */

int iters;							// counter for simulation replications 

double sum_homo;						// sum of the frequencies of the homozygotes 

double rel_major_homo;					// frequency of the major homozygotes among homozygotes 

double poisson[60];					// Poisson probabilities of various coverages 
double denom, sump;					// calculations for the coverage probabilities
double cumpoi[60];					// cumulative Poisson distribution

int maxcoverage;					// maximum coverage level 

double mlD_A;					// candidate disequilibrium coefficient 
double best_error13;				// one third of the estimated error rate

int ig, jg, kg, lg, mg;		// indicators for loops 
int mag, mdg, mgg, mig;
int max_mdg;

int ind_read[nsample+1][5];						// number of observed nucleotides of the four types at a site of an individual 
int pop_read[niters+1][5];						// number of observed nucleotides of the four types at a site in the population sample 
int pop_coverage[niters+1];						// sum of the number of observed nucleotides of the four types at a site in the population sample 

double rancov;						// random number to determine coverage of individual 

int indcov[nsample+1];							// coverage for the individual 

double test_geno;						// random number draw for determining the genotype of the individual 
int geno;							// genotype of the individual (1:major homo, 2:hetero, 3:minor homo) 

double test_nuc;						// random number draw for determining the nucleotide of the erroneous read 

int num_error;						// number of errors at the site 

int num_major_read, num_minor_read;			// number of sequences from major and minor alleles in heterozygotes 

double maxll;						// maximum log likelihood for best estimates 
double llhood;				// log likelihood 
double null_maxll;			// maximum log likelihood under the hypothesis of monomorphism 
double HWE_maxll;			// maximum log likelihood under the hypothesis of HWE
double pol_llstat, HWE_llstat;

int reps_correct_alleles;		// number of simulation replications with correct allele identification 
int s_reps_pol;				// number of simulation replications with MAF greater than zero 
int reps_SNP;				// number of simulation replications with significant pol_llstat at the 5% level 
int reps_nonHWE;			// number of simulation replications with significant HWE_llstat at the 5% level
int reps_pol;				// number of simulation replications with MAF estimates greater than zero 

int n1, n2;							// candidate alleles at the site

int ml_pop_read[niters+1][5];			// number of observed nucleotides of the four types at a site in the population sample for use in the ML estimation 

int ml_ind_read[nsample+1][4];                  // number of observed nucleotides of the three types at a site in an individual for use in the ML estimation 
int null_ind_read[nsample+1][3];			// number of observed nucleotides of the two types at a site of an individual for use in the likelihood function under the null hypothesis 

double null_prob_obs_nuc[nsample+1], HWE_prob_obs_nuc[nsample+1];		// probabilities for use in the likelihood function

double size_grid_D_A, mlD_Amin, mlD_Amax;	// Variables for use in the grid search of the ML estimation of the disequilibrium coefficient D_A

double pre_best_P, pre_best_Q, mlP, mlQ;

double prob_geno[4], prob_nuc[4][4], prob_obs_nuc[nsample+1];		// probabilities for use in the likelihood function 

double best_error, null_best_error, num_best_p, den_best_p, pre_best_p, best_p, best_q, best_D_A, best_f, best_P, best_H, best_Q, best_major_allele_freq, best_minor_allele_freq, best_freq_major_homo, best_freq_minor_homo;			// ML estimates

double freq_major_homo, freq_minor_homo, freq_hetero, error_rate, major_allele_freq, sD_A;

double sum_freq_major_homo, sum_freq_minor_homo, sum_freq_hetero, sum_error_rate, sum_major_allele_freq, sum_sD_A, sum_pol_sD_A, sum_sq_freq_major_homo, sum_sq_freq_minor_homo, sum_sq_freq_hetero, sum_sq_error_rate, sum_sq_major_allele_freq, sum_sq_sD_A, sum_sq_pol_sD_A;
double mean_freq_major_homo, mean_freq_minor_homo, mean_freq_hetero, mean_error_rate, mean_major_allele_freq, mean_sD_A, mean_pol_sD_A;
double mean_sq_freq_major_homo, mean_sq_freq_minor_homo, mean_sq_freq_hetero, mean_sq_error_rate, mean_sq_major_allele_freq, mean_sq_sD_A, mean_sq_pol_sD_A;
double var_freq_major_homo, var_freq_minor_homo, var_freq_hetero, var_error_rate, var_major_allele_freq, var_sD_A, var_pol_sD_A, sd_freq_major_homo, sd_freq_minor_homo, sd_freq_hetero, sd_error_rate, sd_major_allele_freq, sd_sD_A, sd_pol_sD_A;
double sum_inb_coef, sum_pol_inb_coef, sum_sq_inb_coef, sum_sq_pol_inb_coef, inb_coef, mean_inb_coef, mean_sq_inb_coef, var_inb_coef, sd_inb_coef, mean_pol_inb_coef, mean_sq_pol_inb_coef, var_pol_inb_coef, sd_pol_inb_coef;  

double sum_best_freq_major_homo, sum_best_freq_minor_homo, sum_best_H, sum_best_D_A, sum_pol_best_D_A, sum_SNP_best_D_A, sum_best_error, sum_best_major_allele_freq, sum_best_minor_allele_freq, sum_sq_best_freq_major_homo, sum_sq_best_freq_minor_homo, sum_sq_best_H, sum_sq_best_D_A, sum_sq_pol_best_D_A, sum_sq_SNP_best_D_A, sum_sq_best_error, sum_sq_best_major_allele_freq, sum_sq_best_minor_allele_freq;
double mean_best_freq_major_homo, mean_best_freq_minor_homo, mean_best_H, mean_best_D_A, mean_pol_best_D_A, mean_SNP_best_D_A, mean_best_error, mean_best_major_allele_freq, mean_best_minor_allele_freq, msq_best_freq_major_homo, msq_best_freq_minor_homo, msq_best_H, msq_best_D_A, msq_pol_best_D_A, msq_SNP_best_D_A, msq_best_error, msq_best_major_allele_freq, msq_best_minor_allele_freq; 
double var_best_freq_major_homo, var_best_freq_minor_homo, var_best_H, var_best_D_A, var_pol_best_D_A, var_SNP_best_D_A, var_best_error, var_best_major_allele_freq, var_best_minor_allele_freq, sd_best_freq_major_homo, sd_best_freq_minor_homo, sd_best_H, sd_best_D_A, sd_pol_best_D_A, sd_SNP_best_D_A, sd_best_error, sd_best_major_allele_freq, sd_best_minor_allele_freq;
double sum_SD_pol_best_D_A, sum_SD_SNP_best_D_A, mean_SD_pol_best_D_A, mean_SD_SNP_best_D_A, RMSD_pol_best_D_A, RMSD_SNP_best_D_A;
double sum_pol_best_f, sum_sq_pol_best_f, sum_SD_pol_best_f, sum_SNP_best_f, sum_sq_SNP_best_f, sum_SD_SNP_best_f, mean_best_f, sum_best_f, msq_best_f, sum_sq_best_f, var_best_f, sd_best_f, mean_pol_best_f, msq_pol_best_f, var_pol_best_f, sd_pol_best_f, mean_SD_pol_best_f, RMSD_pol_best_f, mean_SNP_best_f, msq_SNP_best_f, var_SNP_best_f, sd_SNP_best_f, mean_SD_SNP_best_f, RMSD_SNP_best_f;

int tot_read, tot_major_homo, tot_minor_homo, tot_hetero, tot_error;

int adjust;	// one if the adjustments are needed, and zero otherwise
double test_best_p, pre_mlD_Amax;
int factor;

// Open the output file. 
stream = fopen( out_file_name, "w" );
if (stream == NULL ) {
	fprintf(stderr, "Cannot open %s for writing.\n", out_file_name); 
	exit(1);
}

// Calculate the Poisson probabilities for various coverages. 

poisson[0] = exp(-double(coverage));
poisson[1] = double(coverage) * exp(-double(coverage));
sump = poisson[0] + poisson[1];

denom = 1.0;

for (ig = 2; ig <= 50; ++ig) { 
	denom = denom * double(ig);
	poisson[ig] = exp(-double(coverage)) * pow(double(coverage),double(ig)) / denom;
	sump = sump + poisson[ig];

	if (sump > 0.999) {
		maxcoverage = ig;
		poisson[ig] = poisson[ig] + (1.0 - sump);
		ig = 1000;
	}
}

printf("\n"); 
printf("Maximum coverage = %d \n", maxcoverage); 
printf("\n"); 
printf("\n"); 

// calculate the cumulative probabilities for coverages
cumpoi[0] = poisson[0];			
for (ig = 1; ig <= maxcoverage; ++ig) {
	cumpoi[ig] = cumpoi[ig-1] + poisson[ig];
}

// Print out the parameter values 
printf("parameters:\n");
fprintf(stream, "parameters:\n");
printf("P = %f\tQ = %f\terror rate = %f\tD_A = %f\taverage coverage per site = %d\tsample size = %d\t# of iterations = %d\n", P, Q, D_A, epsilon, coverage, nsample, niters);
fprintf(stream, "P = %f\tQ = %f\tD_A = %f\terror rate = %f\taverage coverage per site = %d\tsample size = %d\t# of iterations = %d\n", P, Q, D_A, epsilon, coverage, nsample, niters);
printf("\n");
fprintf(stream, "\n");

// Sum of the frequencies of homozygotes, and relative frequency of the major homozygote. 
sum_homo = P + Q;											
rel_major_homo = (double)P/sum_homo;		

// Set the starting values for various quantities.
sum_freq_major_homo = 0.0;
sum_freq_minor_homo = 0.0;
sum_freq_hetero = 0.0;
sum_major_allele_freq = 0.0;
sum_error_rate = 0.0;
sum_inb_coef = 0.0;
sum_pol_inb_coef = 0.0;
sum_sD_A = 0.0;
sum_pol_sD_A = 0.0;
sum_sq_freq_major_homo = 0.0;
sum_sq_freq_minor_homo = 0.0;
sum_sq_freq_hetero = 0.0;
sum_sq_major_allele_freq = 0.0;
sum_sq_error_rate = 0.0;
sum_sq_inb_coef = 0.0;
sum_sq_pol_inb_coef = 0.0;
sum_sq_sD_A = 0.0;
sum_sq_pol_sD_A = 0.0;

reps_correct_alleles = 0;
s_reps_pol = 0;
reps_pol = 0;
reps_SNP = 0;
reps_nonHWE = 0;
sum_best_freq_major_homo = 0.0;
sum_best_major_allele_freq = 0.0;
sum_best_H = 0.0;
sum_best_freq_minor_homo = 0.0;
sum_best_minor_allele_freq = 0.0;
sum_best_D_A = 0.0;
sum_pol_best_D_A = 0.0;
sum_SNP_best_D_A = 0.0;
sum_best_f = 0.0;
sum_pol_best_f = 0.0;
sum_SNP_best_f = 0.0;					
sum_best_error = 0.0;
sum_sq_best_freq_major_homo = 0.0;
sum_sq_best_H = 0.0;
sum_sq_best_freq_minor_homo = 0.0;
sum_sq_best_major_allele_freq = 0.0;
sum_sq_best_minor_allele_freq = 0.0;
sum_sq_best_D_A = 0.0;
sum_sq_pol_best_D_A = 0.0;
sum_sq_SNP_best_D_A = 0.0;
sum_sq_best_f = 0.0;
sum_sq_pol_best_f = 0.0;
sum_sq_SNP_best_f = 0.0;
sum_sq_best_error = 0.0;
sum_SD_pol_best_D_A = 0.0;
sum_SD_SNP_best_D_A = 0.0;
sum_SD_pol_best_f = 0.0;
sum_SD_SNP_best_f = 0.0;


/* *********************************************************************************************************************** */

// Print out the field names
fprintf(stream, "iteration\tsample P\tsample Q\tsample error rate\tsample p\tsample D_A\tbest_major_allele_freq\tbest_D_A\tbest_f\testimated P\testimated Q\tbest_error\tHWE_maxll\tnull_maxll\tmaxll\tpol_llstat\tHWE_llstat\n");
printf("iteration\tsample P\tsample Q\tsample error rate\tsample p\tsample D_A\tbest_major_allele_freq\tbest_D_A\tbest_f\testimated P\testimated Q\tbest_error\tHWE_maxll\tnull_maxll\tmaxll\tpol_llstat\tHWE_llstat\n");

// Loop over simulations for different population samples, in each case generating estimates of genotyope frequencies.

for (iters = 1; iters <= niters; ++iters) {

	tot_read = 0;
	tot_major_homo = 0;
	tot_minor_homo = 0;
	tot_hetero = 0;
	tot_error = 0;


	// Zero the counters for the sequence reads in the sample.
	for (ig = 1; ig <= 4; ++ig){  
		pop_read[iters][ig] = 0;
	}

	// Start generating the read data for each individual.
	for (jg = 1; jg <= nsample; ++jg) {		

		rancov = ranmar();								// Determine coverage for the individual at the site.
		for (ig = 0; ig <= maxcoverage; ++ig) {
			if (rancov <= cumpoi[ig]) {
				indcov[jg] = ig;
				ig = 10000;
			}
		}

		tot_read = tot_read + indcov[jg];

		for (mg = 1; mg <= 4; ++mg)							// Initialize the read counts of the individual to zero 
			ind_read[jg][mg] = 0;							

		if (ranmar() < sum_homo) {						// Generate sequence-read data when homozygous 

			test_geno = ranmar();							
			if (test_geno < rel_major_homo) {					// Major homozygote 
				geno = 1;
				tot_major_homo = tot_major_homo + 1;
				num_error = ignbin(indcov[jg],epsilon);                             // number of erroneous reads for this individual 
				tot_error = tot_error + num_error;
				ind_read[jg][1] = ind_read[jg][1] + indcov[jg] - num_error;                       // number of correct reads for this individual 
				if (num_error > 0) {
                                	for (kg = 1; kg <= num_error; ++kg) {
                                        	test_nuc = ranmar();
                                        	if (test_nuc < 0.3333333){
                                                	ind_read[jg][2] = ind_read[jg][2] + 1;
                                        	} else if (test_nuc < 0.6666667) {
                                                	ind_read[jg][3] = ind_read[jg][3] + 1;
                                        	} else {
                                                	ind_read[jg][4] = ind_read[jg][4] + 1;
                                        	}
                                	}
                        	}
			} else {								// Minor homozygote 
				geno = 3;
				tot_minor_homo = tot_minor_homo + 1;
				num_error = ignbin(indcov[jg],epsilon);                             // number of erroneous reads for this individual 
				tot_error = tot_error + num_error;
				ind_read[jg][2] = ind_read[jg][2] + indcov[jg] - num_error;                       // number of correct reads for this individual 
				if (num_error > 0) {
                                	for (kg = 1; kg <= num_error; ++kg) {
                                        	test_nuc = ranmar();
                                        	if (test_nuc < 0.3333333){
                                                	ind_read[jg][1] = ind_read[jg][1] + 1;
                                        	} else if (test_nuc < 0.6666667) {
                                                	ind_read[jg][3] = ind_read[jg][3] + 1;
                                        	} else {
                                                	ind_read[jg][4] = ind_read[jg][4] + 1;
                                        	}
                                	}
                        	}
			}


		} else {					// generate sequence read data when heterozygous 
			geno = 2;
			num_major_read = ignbin(indcov[jg],0.5);					// number of draws of the major allele 
			num_minor_read = indcov[jg] - num_major_read;					// number of draws of the minor allele 
			tot_hetero = tot_hetero + 1;
			
			if (num_major_read >= 1) {
				num_error = ignbin(num_major_read,epsilon);
				tot_error = tot_error + num_error;
				ind_read[jg][1] = ind_read[jg][1] + num_major_read - num_error;
				if (num_error > 0) {
					for (kg = 1; kg <= num_error; ++kg) {
						test_nuc = ranmar();
						if (test_nuc < 0.3333333){
							ind_read[jg][2] = ind_read[jg][2] + 1;
						} else if (test_nuc < 0.6666667) {
							ind_read[jg][3] = ind_read[jg][3] + 1;
						} else {
							ind_read[jg][4] = ind_read[jg][4] + 1;
						}
					}
				}
			}

			if (num_minor_read >= 1) {
				num_error = ignbin(num_minor_read,epsilon);
				tot_error = tot_error + num_error;
				ind_read[jg][2] = ind_read[jg][2] + num_minor_read - num_error;
				if (num_error > 0) {
					for (kg = 1; kg <= num_error; ++kg) {
						test_nuc = ranmar();
						if (test_nuc < 0.3333333) {
							ind_read[jg][1] = ind_read[jg][1] + 1;
						} else if (test_nuc < 0.6666667) {
							ind_read[jg][3] = ind_read[jg][3] + 1;
						} else {
							ind_read[jg][4] = ind_read[jg][4] + 1;
						}
					}
				}
			}
		}
                for(lg=1; lg<=4; ++lg){
                	pop_read[iters][lg] = pop_read[iters][lg] + ind_read[jg][lg];
		}
	}		// end of the loop over the individuals
	freq_major_homo = (double(tot_major_homo)/double(nsample));
	freq_minor_homo = (double(tot_minor_homo)/double(nsample));
	freq_hetero = (double(tot_hetero)/double(nsample));
	error_rate = (double(tot_error)/double(tot_read));
	major_allele_freq = freq_major_homo + freq_hetero/2.0;
	if (major_allele_freq < 1.0) {
		s_reps_pol = s_reps_pol + 1;
		sD_A = major_allele_freq*(1.0-major_allele_freq) - freq_hetero/2.0;
		sum_pol_sD_A = sum_pol_sD_A + sD_A;
		sum_sq_pol_sD_A = sum_sq_pol_sD_A + sD_A*sD_A;
		inb_coef = 1.0 - freq_hetero/( 2.0*major_allele_freq*(1.0-major_allele_freq) );
		sum_pol_inb_coef = sum_pol_inb_coef + inb_coef;
		sum_sq_pol_inb_coef = sum_sq_pol_inb_coef + inb_coef*inb_coef;
	} else {
		sD_A = 0.0;
		inb_coef = 0.0;
	}
	sum_freq_major_homo = sum_freq_major_homo + freq_major_homo;
	sum_freq_minor_homo = sum_freq_minor_homo + freq_minor_homo;
	sum_freq_hetero = sum_freq_hetero + freq_hetero;
	sum_error_rate = sum_error_rate + error_rate;
	sum_major_allele_freq = sum_major_allele_freq + major_allele_freq;
	sum_sD_A = sum_sD_A + sD_A;
	sum_inb_coef = sum_inb_coef + inb_coef;
	sum_sq_freq_major_homo = sum_sq_freq_major_homo + freq_major_homo*freq_major_homo;
        sum_sq_freq_minor_homo = sum_sq_freq_minor_homo + freq_minor_homo*freq_minor_homo;
        sum_sq_freq_hetero = sum_sq_freq_hetero + freq_hetero*freq_hetero;
        sum_sq_error_rate = sum_sq_error_rate + error_rate*error_rate;
	sum_sq_major_allele_freq = sum_sq_major_allele_freq + major_allele_freq*major_allele_freq;
	sum_sq_sD_A = sum_sq_sD_A + sD_A*sD_A;
	sum_sq_inb_coef = sum_sq_inb_coef + inb_coef*inb_coef;
	
/* ******************************************************************************************************************** */

	// Obtain the ML estimates of the frequencies of homozygotes and error rate for this sample.
	maxll = -10000000000.0;

	// Find the most and second most abundant nucleotides at the site and consider them as potential alleles at the site.
	if (pop_read[iters][1] >= pop_read[iters][2]) {
		n1 = 1;
		n2 = 2;
	} else {
		n1 = 2;
		n2 = 1;
	}
	if (pop_read[iters][3] >= pop_read[iters][n1]) {
		n2 = n1;
		n1 = 3;
	} else if (pop_read[iters][3] >= pop_read[iters][n2]) {
		n2 = 3;
	}
	if (pop_read[iters][4] >= pop_read[iters][n1]) {
		n2 = n1;
		n1 = 4;
	} else if (pop_read[iters][4] >= pop_read[iters][n2]) {
		n2 = 4;
	}

	pop_coverage[iters] = pop_read[iters][1] + pop_read[iters][2] + pop_read[iters][3] + pop_read[iters][4];

	if (pop_read[iters][n1] < pop_coverage[iters]) {
		if (n1*n2 == 2) {
                	reps_correct_alleles = reps_correct_alleles + 1;
	        }
		/* Count the number of different types of reads at the site for use in the ML estimation */
		ml_pop_read[iters][1] = pop_read[iters][n1];
		ml_pop_read[iters][2] = pop_read[iters][n2];
		ml_pop_read[iters][3] = pop_coverage[iters] - ml_pop_read[iters][1] - ml_pop_read[iters][2];
		for(mig = 1; mig <= nsample; ++mig) {
			ml_ind_read[mig][1] = ind_read[mig][n1];
			ml_ind_read[mig][2] = ind_read[mig][n2];
			ml_ind_read[mig][3] = indcov[mig] - ml_ind_read[mig][1] - ml_ind_read[mig][2];
			null_ind_read[mig][1] = ind_read[mig][n1];
			null_ind_read[mig][2] = indcov[mig] - null_ind_read[mig][1];
		}
		best_error = 1.5*((double)ml_pop_read[iters][3]/(double)pop_coverage[iters]);
		num_best_p = 2.0*ml_pop_read[iters][1] - ml_pop_read[iters][3];
		den_best_p = 2.0*(ml_pop_read[iters][1]+ml_pop_read[iters][2]-ml_pop_read[iters][3]);
		pre_best_p = num_best_p/den_best_p;
		if ( 2.0*nsample*pre_best_p - (int)(2.0*nsample*pre_best_p) >= 0.5 ) {
			best_p = (double)(1.0/(2.0*nsample))*( (int)(2.0*nsample*pre_best_p)+1 );
		} else {
        		best_p = (double)(1.0/(2.0*nsample))*( (int)(2.0*nsample*pre_best_p) );
		}
		if (best_p > 1.0 - 1.0e-8) {
			best_p = (double)(1.0/(2.0*nsample))*( (int)(2.0*nsample*pre_best_p) );
		}
		// Calculate the corresponding maximum log-likelihood under the null hypothesis of monomorphism 
		null_best_error = (pop_coverage[iters]-ml_pop_read[iters][1])/(double)pop_coverage[iters];
		// Sum the log-likelihoods over the individuals
		null_maxll = 0.0;
		for(mig = 1; mig <= nsample; ++mig) {
			if (indcov[mig] > 0) {
				null_prob_obs_nuc[mig] = pow(1.0-null_best_error,(double)null_ind_read[mig][1])*pow(null_best_error/3.0,(double)null_ind_read[mig][2]);
				if (null_prob_obs_nuc[mig] > 0.0) {
					null_maxll = null_maxll + log(null_prob_obs_nuc[mig]);
				} else {
					null_maxll = -10000000001.0;
				}
			}
		}
		
		// Grid search for the preliminary ML estimate of D_A starts here.
		size_grid_D_A = 1.0/(double)nsample;
		if ( -pow(best_p,2.0) >= -pow(1.0-best_p,2.0) ) {
			mlD_Amin = -pow(best_p,2.0);
		} else {
			mlD_Amin = -pow(1.0-best_p,2.0);
		}
		test_best_p = best_p*(double)nsample-(int)(best_p*nsample);
		if (test_best_p < 1.0e-8) {
			mlD_Amax = best_p*(1.0-best_p);
		} else {
			pre_mlD_Amax = best_p*(1.0-best_p);
			factor = (int)( (pre_mlD_Amax-mlD_Amin)/size_grid_D_A );
			mlD_Amax = mlD_Amin + size_grid_D_A*factor;
		}
		mlD_A = mlD_Amin - size_grid_D_A;
		max_mdg = (int)( nsample*(mlD_Amax-mlD_Amin) ) + 2;
		// Loop over the candidate disequilibrium coefficients at the site 
		for (mdg = 1; mdg <= max_mdg; ++mdg) {
			if (mdg == max_mdg) {
				mlD_A = mlD_Amax;
			} else {
				mlD_A = mlD_A + size_grid_D_A;
			}		
			// Calculate the probabilities for use in the likelihood function 
			best_error13 = best_error/3.0;
			// Calculate the probability of each of the three genotypes for an individual
			prob_geno[1] = pow(best_p,2.0) + mlD_A;	// frequency of the homozygotes for the most abundant nucleotide at the site 
			prob_geno[2] = 2.0*( best_p*(1.0-best_p)-mlD_A );      // frequency of heterozygotes at the site 
			prob_geno[3] = pow(1.0-best_p,2.0) + mlD_A;      // frequency of the homozygotes for the second most abundant nucleotide at the site 
			// Calculate the probability of a nucleotide given a genotype of the individual at the site 
			prob_nuc[1][1] = 1.0 - best_error;
			prob_nuc[1][2] = best_error13;
			prob_nuc[1][3] = best_error13;
			prob_nuc[2][1] = (1/2.0)*(1.0-best_error + best_error13);
			prob_nuc[2][2] = (1/2.0)*(best_error13 + 1.0-best_error);
			prob_nuc[2][3] = best_error13;
			prob_nuc[3][1] = best_error13;
			prob_nuc[3][2] = 1.0 - best_error;
			prob_nuc[3][3] = best_error13;
				
			// Sum the log-likelihoods over the individuals 
			llhood = 0.0;
			for(mig = 1; mig <= nsample; ++mig) {
				if (indcov[mig] > 0) {
					// Sum the probabilities over the genotypes of the individual 
					prob_obs_nuc[mig] = 0.0;
					for(mgg = 1; mgg <= 3; ++mgg) {
						prob_obs_nuc[mig] = prob_obs_nuc[mig] + prob_geno[mgg]*pow(prob_nuc[mgg][1],(double)ml_ind_read[mig][1])*pow(prob_nuc[mgg][2],(double)ml_ind_read[mig][2])*pow(prob_nuc[mgg][3],(double)ml_ind_read[mig][3]);
					}
					if (prob_obs_nuc[mig] > 0.0) {
						llhood = llhood + log(prob_obs_nuc[mig]);
					} else {
						llhood = -10000000001.0;
					}
				}
			}

			// Examine whether this is a new ML solution for the sample 
			if (llhood > maxll) {
                        	maxll = llhood;
				best_D_A = mlD_A;
                	}
		}	// End of the loop over the candidate disequilibrium coefficients at the site
		pre_best_P = pow(best_p,2.0) + best_D_A;
		pre_best_Q = pow(1.0-best_p,2.0) + best_D_A;
		if (pre_best_P < 0.0) {
			best_P = 0.0;
		} else {
			if ( nsample*pre_best_P - (int)(nsample*pre_best_P) >= 0.5 ) {
				best_P = (double)(1.0/(double)nsample)*( (int)(nsample*pre_best_P)+1 );
			} else {
				best_P = (double)(1.0/(double)nsample)*( (int)(nsample*pre_best_P) );
			}
		}
		if (pre_best_Q < 0.0) {
			best_Q = 0.0;
		} else {
			if ( nsample*pre_best_Q - (int)(nsample*pre_best_Q) >= 0.5 ) {
				best_Q = (double)(1.0/(double)nsample)*( (int)(nsample*pre_best_Q)+1 );
			} else {
				best_Q = (double)(1.0/(double)nsample)*( (int)(nsample*pre_best_Q) );
			}
		}	

		// The adjustments of the ML genotype-frequency estimates start here.	
		// Adjust the ML estimates by a grid search
		adjust = 1;
		while (adjust == 1) {
			adjust = 0;
			for (mag=1; mag<=8; ++mag) {
				if (mag == 1) {
					mlP = best_P - size_grid_D_A;
					mlQ = best_Q - size_grid_D_A;
				} else if (mag == 2) {
					mlP = best_P - size_grid_D_A;
					mlQ = best_Q;
				} else if (mag == 3) {
					mlP = best_P - size_grid_D_A;
					mlQ = best_Q + size_grid_D_A;
				} else if (mag == 4) {
					mlP = best_P;
					mlQ = best_Q - size_grid_D_A;
				} else if (mag == 5) {
					mlP = best_P;
					mlQ = best_Q + size_grid_D_A;
				} else if (mag == 6) {
					mlP = best_P + size_grid_D_A;
					mlQ = best_Q - size_grid_D_A;
				} else if (mag == 7) {
					mlP = best_P + size_grid_D_A;
					mlQ = best_Q;
				} else if (mag == 8) {
					mlP = best_P + size_grid_D_A;
					mlQ = best_Q + size_grid_D_A;
				}
				if (mlP < -1.0e-8 || mlP > 1.0 + 1.0e-8 || mlQ < -1.0e-8 || mlQ > 1.0 + 1.0e-8 || mlP + mlQ > 1.0 + 1.0e-8) {
					mlP = best_P;
					mlQ = best_Q;
					llhood = -10000000001.0;	
				} else {
					llhood = 0.0;
					prob_geno[1] = mlP;
					prob_geno[2] = (double)1.0 - mlP - mlQ;
					prob_geno[3] = mlQ;
					// Sum the log-likelihoods over the individuals
					for (mig = 1; mig <= nsample; mig++) {
						if (indcov[mig] > 0) {
							// Sum the probabilities over the genotypes of the individual
							prob_obs_nuc[mig] = 0.0;
							for (mgg = 1; mgg <= 3; ++mgg) {
								prob_obs_nuc[mig] = prob_obs_nuc[mig] + prob_geno[mgg]*pow(prob_nuc[mgg][1], (double)ml_ind_read[mig][1])*pow(prob_nuc[mgg][2], (double)ml_ind_read[mig][2])*pow(prob_nuc[mgg][3], (double)ml_ind_read[mig][3]);
							}
							if (prob_obs_nuc[mig] > 0.0) {
								llhood = llhood + log(prob_obs_nuc[mig]);
							} else {
								llhood = -10000000001.0;
							}
						}
					}
					if (llhood > maxll) {
						maxll = llhood;
						adjust = 1;
						best_P = mlP;
						best_Q = mlQ;
					}
				} 
			}
		}
		best_p = best_P + ( (double)1.0-best_P-best_Q )/(double)2.0;
		// Calculate the maximum log-likelihood under the hypothesis of HWE
		HWE_maxll = 0.0;
		for (mig = 1; mig <= nsample; mig++) {
			if (indcov[mig] > 0) {
				HWE_prob_obs_nuc[mig] = pow( best_p, 2.0 )*pow( 1.0-best_error, (double)ml_ind_read[mig][1] )*pow( best_error/3.0, (double)ml_ind_read[mig][2] )*pow( best_error/3.0, (double)ml_ind_read[mig][3] ) + 2.0*best_p*(1.0-best_p)*pow( 0.5-best_error/3.0, (double)(ml_ind_read[mig][1]+ml_ind_read[mig][2]) )*pow( best_error/3.0, (double)ml_ind_read[mig][3] ) + pow( 1.0-best_p, 2.0 )*pow( best_error/3.0, (double)ml_ind_read[mig][1] )*pow( 1.0-best_error, (double)ml_ind_read[mig][2] )*pow( best_error/3.0, (double)ml_ind_read[mig][3] );
				if (HWE_prob_obs_nuc[mig] > 0.0) {
					HWE_maxll = HWE_maxll + log(HWE_prob_obs_nuc[mig]);
				} else {
					HWE_maxll = -10000000001.0;
				}
			}
		}
		if (HWE_maxll >= maxll) {
			maxll = HWE_maxll;
		}
		if (null_maxll >= maxll) {
			maxll = null_maxll;
			HWE_maxll = null_maxll;
			best_P = 1.0;
			best_Q = 0.0;
			best_error = null_best_error;
		}
		best_H = (double)1.0 - best_P - best_Q;
		best_p = best_P + best_H/(double)2.0;				
		best_q = 1.0 - best_p;
		best_D_A = best_P - pow(best_p,2.0);
		if ( best_p <= 1.0-(double)1.0/(2.0*nsample) ) {
			best_f = (2.0*best_p*best_q-best_H)/(2.0*best_p*best_q);
		} else {
			best_f = 0;
		}
		if (n1 == 1) {
			best_major_allele_freq = best_p;
			best_minor_allele_freq = best_q;
        		best_freq_major_homo = best_P;
                	best_freq_minor_homo = best_Q;
        	} else if (n1 == 2) {
			best_major_allele_freq = best_q;
			best_minor_allele_freq = best_p;
                	best_freq_major_homo = best_Q;
                	best_freq_minor_homo = best_P;
        	}
		if (best_p < 1.0) {
			reps_pol = reps_pol + 1;
			sum_pol_best_D_A = sum_pol_best_D_A + best_D_A;
			sum_sq_pol_best_D_A = sum_sq_pol_best_D_A + best_D_A*best_D_A;
			sum_SD_pol_best_D_A = sum_SD_pol_best_D_A + pow(best_D_A - D_A, 2.0);
			sum_pol_best_f = sum_pol_best_f + best_f;
                        sum_sq_pol_best_f = sum_sq_pol_best_f + best_f*best_f;
                        sum_SD_pol_best_f = sum_SD_pol_best_f + pow(best_f - f, 2.0);
		}
		pol_llstat = 2.0*(maxll - null_maxll);
		HWE_llstat = 2.0*(maxll - HWE_maxll);
		if (pol_llstat > 5.991) {
			reps_SNP = reps_SNP + 1;
			sum_SNP_best_D_A = sum_SNP_best_D_A + best_D_A;
			sum_sq_SNP_best_D_A = sum_sq_SNP_best_D_A + best_D_A*best_D_A;
			sum_SD_SNP_best_D_A = sum_SD_SNP_best_D_A + pow(best_D_A - D_A, 2.0);
			sum_SNP_best_f = sum_SNP_best_f + best_f;
                        sum_sq_SNP_best_f = sum_sq_SNP_best_f + best_f*best_f;
                        sum_SD_SNP_best_f = sum_SD_SNP_best_f + pow(best_f - f, 2.0);
			if (HWE_llstat > 3.841) {
                        	reps_nonHWE = reps_nonHWE + 1;
                	}
		}			 
		fprintf(stream, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", iters, freq_major_homo, freq_minor_homo, error_rate, major_allele_freq, sD_A, best_major_allele_freq, best_D_A, best_f, best_freq_major_homo, best_freq_minor_homo, best_error, HWE_maxll, null_maxll, maxll, pol_llstat, HWE_llstat);
		printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", iters, freq_major_homo, freq_minor_homo, error_rate, major_allele_freq, sD_A, best_major_allele_freq, best_D_A, best_f, best_freq_major_homo, best_freq_minor_homo, best_error, HWE_maxll, null_maxll, maxll, pol_llstat, HWE_llstat);
	} else if (pop_read[iters][n1] == pop_coverage[iters]) {
		if (n1 == 1) {
			reps_correct_alleles = reps_correct_alleles + 1;
		}
		best_freq_major_homo = 1.0;
		best_freq_minor_homo = 0.0;
		best_H = 0.0;
		best_p = 1.0;
		best_q = 0.0;
		best_error = 0.0;
		best_D_A = 0.0;
		best_f = 0.0;
		best_major_allele_freq = best_p;
		best_minor_allele_freq = best_q;
		pol_llstat = 0.0;
		HWE_llstat = 0.0;
		fprintf(stream, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\tNA\tNA\tNA\t%f\t%f\n", iters, freq_major_homo, freq_minor_homo, error_rate, major_allele_freq, sD_A, best_major_allele_freq, best_D_A, best_f, best_freq_major_homo, best_freq_minor_homo, best_error, pol_llstat, HWE_llstat);
		printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\tNA\tNA\tNA\t%f\t%f\n", iters, freq_major_homo, freq_minor_homo, error_rate, major_allele_freq, sD_A, best_major_allele_freq, best_D_A, best_f, best_freq_major_homo, best_freq_minor_homo, best_error, pol_llstat, HWE_llstat);     
	}		    
	sum_best_freq_major_homo = sum_best_freq_major_homo + best_freq_major_homo;					// increment the summary statistics 
	sum_best_freq_minor_homo = sum_best_freq_minor_homo + best_freq_minor_homo;
	sum_best_H = sum_best_H + best_H;
	sum_best_D_A = sum_best_D_A + best_D_A;
	sum_best_f = sum_best_f + best_f;
	
	sum_best_error = sum_best_error + best_error;
	sum_best_major_allele_freq = sum_best_major_allele_freq + best_major_allele_freq;
	sum_best_minor_allele_freq = sum_best_minor_allele_freq + best_minor_allele_freq;
	sum_sq_best_freq_major_homo = sum_sq_best_freq_major_homo + pow(best_freq_major_homo,2.0);
	sum_sq_best_freq_minor_homo = sum_sq_best_freq_minor_homo + pow(best_freq_minor_homo,2.0);
	sum_sq_best_H = sum_sq_best_H + pow(best_H,2.0);
	sum_sq_best_D_A = sum_sq_best_D_A + pow(best_D_A,2.0);
	sum_sq_best_f = sum_sq_best_f + pow(best_f,2.0);
	sum_sq_best_error = sum_sq_best_error + pow(best_error,2.0);
	sum_sq_best_major_allele_freq = sum_sq_best_major_allele_freq + pow(best_major_allele_freq,2.0);
	sum_sq_best_minor_allele_freq = sum_sq_best_minor_allele_freq + pow(best_minor_allele_freq,2.0); 
}			// end of the loop for generating the reads over the population samples 

// Calculate the means and variances of estimates over simulation replications
mean_best_freq_major_homo = sum_best_freq_major_homo / double(niters);
mean_best_freq_minor_homo = sum_best_freq_minor_homo / double(niters);
mean_best_H = sum_best_H / double(niters);
mean_best_D_A = sum_best_D_A / double(niters);
mean_best_f = sum_best_f / double(niters);
mean_best_error = sum_best_error / double(niters);
mean_best_major_allele_freq = sum_best_major_allele_freq / double(niters);
mean_best_minor_allele_freq = sum_best_minor_allele_freq / double(niters);
msq_best_freq_major_homo = sum_sq_best_freq_major_homo / double(niters);
msq_best_freq_minor_homo = sum_sq_best_freq_minor_homo / double(niters);
msq_best_H = sum_sq_best_H / double(niters);
msq_best_D_A = sum_sq_best_D_A / double(niters);
msq_best_f = sum_sq_best_f / double(niters);
msq_best_error = sum_sq_best_error / double(niters);
msq_best_major_allele_freq = sum_sq_best_major_allele_freq / double(niters);
msq_best_minor_allele_freq = sum_sq_best_minor_allele_freq / double(niters);
var_best_freq_major_homo = (double(niters)/(double(niters) - 1.0)) * (msq_best_freq_major_homo - pow(mean_best_freq_major_homo,2.0));
var_best_freq_minor_homo = (double(niters)/(double(niters) - 1.0)) * (msq_best_freq_minor_homo - pow(mean_best_freq_minor_homo,2.0));
var_best_H = (double(niters)/(double(niters) - 1.0)) * (msq_best_H - pow(mean_best_H,2.0));
var_best_D_A = (double(niters)/(double(niters) - 1.0)) * (msq_best_D_A - pow(mean_best_D_A,2.0));
var_best_f = (double(niters)/(double(niters) - 1.0)) * (msq_best_f - pow(mean_best_f,2.0));
var_best_error = (double(niters)/(double(niters) - 1.0)) * (msq_best_error - pow(mean_best_error,2.0));
var_best_major_allele_freq = (double(niters)/(double(niters) - 1.0)) * (msq_best_major_allele_freq - pow(mean_best_major_allele_freq,2.0));
var_best_minor_allele_freq = (double(niters)/(double(niters) - 1.0)) * (msq_best_minor_allele_freq - pow(mean_best_minor_allele_freq,2.0));
sd_best_freq_major_homo = pow(var_best_freq_major_homo,0.5);
sd_best_freq_minor_homo = pow(var_best_freq_minor_homo,0.5);
sd_best_H = pow(var_best_H,0.5);
sd_best_D_A = pow(var_best_D_A,0.5);
sd_best_f = pow(var_best_f,0.5);
sd_best_error = pow(var_best_error,0.5);
sd_best_major_allele_freq = pow(var_best_major_allele_freq,0.5);
sd_best_minor_allele_freq = pow(var_best_minor_allele_freq,0.5);
if (reps_pol > 0) {
	mean_pol_best_D_A = sum_pol_best_D_A / double(reps_pol);
	msq_pol_best_D_A = sum_sq_pol_best_D_A / double(reps_pol);
	var_pol_best_D_A = (double(reps_pol)/(double(reps_pol) - 1.0)) * (msq_pol_best_D_A - pow(mean_pol_best_D_A,2.0));
	sd_pol_best_D_A = pow(var_pol_best_D_A,0.5);
	mean_SD_pol_best_D_A = sum_SD_pol_best_D_A / double(reps_pol);
	RMSD_pol_best_D_A = pow(mean_SD_pol_best_D_A,0.5);
	mean_pol_best_f = sum_pol_best_f / double(reps_pol);
        msq_pol_best_f = sum_sq_pol_best_f / double(reps_pol);
        var_pol_best_f = (double(reps_pol)/(double(reps_pol) - 1.0)) * (msq_pol_best_f - pow(mean_pol_best_f,2.0));
        sd_pol_best_f = pow(var_pol_best_f,0.5);
        mean_SD_pol_best_f = sum_SD_pol_best_f / double(reps_pol);
        RMSD_pol_best_f = pow(mean_SD_pol_best_f,0.5);
}
if (reps_SNP > 0) {
        mean_SNP_best_D_A = sum_SNP_best_D_A / double(reps_SNP);
        msq_SNP_best_D_A = sum_sq_SNP_best_D_A / double(reps_SNP);
        var_SNP_best_D_A = (double(reps_SNP)/(double(reps_SNP) - 1.0)) * (msq_SNP_best_D_A - pow(mean_SNP_best_D_A,2.0));
        sd_SNP_best_D_A = pow(var_SNP_best_D_A,0.5);
	mean_SD_SNP_best_D_A = sum_SD_SNP_best_D_A / double(reps_SNP);
        RMSD_SNP_best_D_A = pow(mean_SD_SNP_best_D_A,0.5);
	mean_SNP_best_f = sum_SNP_best_f / double(reps_SNP);
        msq_SNP_best_f = sum_sq_SNP_best_f / double(reps_SNP);
        var_SNP_best_f = (double(reps_SNP)/(double(reps_SNP) - 1.0)) * (msq_SNP_best_f - pow(mean_SNP_best_f,2.0));
        sd_SNP_best_f = pow(var_SNP_best_f,0.5);
        mean_SD_SNP_best_f = sum_SD_SNP_best_f / double(reps_SNP);
        RMSD_SNP_best_f = pow(mean_SD_SNP_best_f,0.5);
}

fprintf(stream, "\n");
fprintf(stream, "\n"); 
printf("\n");
printf("\n");
fprintf(stream, "Estimation results\n");
printf("Estimation results\n");
fprintf(stream, "Rate of correct allele identification: %f\n", (double)reps_correct_alleles/(double)niters);
printf("Rate of correct allele identification: %f\n", (double)reps_correct_alleles/(double)niters);
fprintf(stream, "Rate of significant polyrmophism detection: %f\n", (double)reps_SNP/(double)niters);
printf("Rate of significant polyrmophism detection: %f\n", (double)reps_SNP/(double)niters);
if (reps_SNP > 0) {
	fprintf(stream, "Rate of HWE-deviation detection: %f\n", (double)reps_nonHWE/(double)reps_SNP);
	printf("Rate of HWE-deviation detection: %f\n", (double)reps_nonHWE/(double)reps_SNP);
} else {
	fprintf(stream, "Rate of HWE-deviation detection: NA\n");
	printf("Rate of HWE-deviation detection: NA\n");
}
fprintf(stream, "quantity\taverage\tSD\tRMSD\tnumber of data\n");
printf("quantity\taverage\tSD\tRMSD\tnumber of data\n");
fprintf(stream, "major allele frequency\t%6.5f\t%6.5f\tNA\t%d\n", mean_best_major_allele_freq, sd_best_major_allele_freq, niters);
printf("major allele frequency\t%6.5f\t%6.5f\tNA\t%d\n", mean_best_major_allele_freq, sd_best_major_allele_freq, niters);
fprintf(stream, "minor allele frequency\t%6.5f\t%6.5f\tNA\t%d\n", mean_best_minor_allele_freq, sd_best_minor_allele_freq, niters);
printf("minor allele frequency\t%6.5f\t%6.5f\tNA\t%d\n", mean_best_minor_allele_freq, sd_best_minor_allele_freq, niters);
fprintf(stream, "disequilibrium coefficient\t%6.5f\t%6.5f\tNA\t%d\n", mean_best_D_A, sd_best_D_A, niters);
printf("disequilibrium coefficient\t%6.5f\t%6.5f\tNA\t%d\n", mean_best_D_A, sd_best_D_A, niters);
fprintf(stream, "inbreeding coefficient\t%6.5f\t%6.5f\tNA\t%d\n", mean_best_f, sd_best_f, niters);
printf("inbreeding coefficient\t%6.5f\t%6.5f\tNA\t%d\n", mean_best_f, sd_best_f, niters);
if (reps_pol > 0) {
        fprintf(stream, "polymorphic disequilibrium coefficient\t%6.5f\t%6.5f\t%6.5f\t%d\n", mean_pol_best_D_A, sd_pol_best_D_A, RMSD_pol_best_D_A, reps_pol);
        printf("polymorphic disequilibrium coefficient\t%6.5f\t%6.5f\t%6.5f\t%d\n", mean_pol_best_D_A, sd_pol_best_D_A, RMSD_pol_best_D_A, reps_pol);
	fprintf(stream, "polymorphic inbreeding coefficient\t%6.5f\t%6.5f\t%6.5f\t%d\n", mean_pol_best_f, sd_pol_best_f, RMSD_pol_best_f, reps_pol);
        printf("polymorphic inbreeding coefficient\t%6.5f\t%6.5f\t%6.5f\t%d\n", mean_pol_best_f, sd_pol_best_f, RMSD_pol_best_f, reps_pol);
} else {
        fprintf(stream, "polymorphic disequiibrium coefficient\tNA\tNA\tNA\t%d\n", reps_pol);
        printf("polymorphic disequilibrium coefficient\tNA\tNA\tNA\t%d\n", reps_pol);
	fprintf(stream, "polymorphic inbreeding coefficient\tNA\tNA\tNA\t%d\n", reps_pol);
        printf("polymorphic inbreeding coefficient\tNA\tNA\tNA\t%d\n", reps_pol);
}
if (reps_SNP > 0) {
        fprintf(stream, "SNP disequilibrium coefficient\t%6.5f\t%6.5f\t%6.5f\t%d\n", mean_SNP_best_D_A, sd_SNP_best_D_A, RMSD_SNP_best_D_A, reps_SNP);
        printf("SNP disequilibrium coefficient\t%6.5f\t%6.5f\t%6.5f\t%d\n", mean_SNP_best_D_A, sd_SNP_best_D_A, RMSD_SNP_best_D_A, reps_SNP);
	fprintf(stream, "SNP inbreeding coefficient\t%6.5f\t%6.5f\t%6.5f\t%d\n", mean_SNP_best_f, sd_SNP_best_f, RMSD_SNP_best_f, reps_SNP);
        printf("SNP inbreeding coefficient\t%6.5f\t%6.5f\t%6.5f\t%d\n", mean_SNP_best_f, sd_SNP_best_f, RMSD_SNP_best_f, reps_SNP);
} else {
        fprintf(stream, "SNP disequilibrium coefficient\tNA\tNA\t%d\n", reps_SNP);
        printf("SNP disequilibrium coefficient\tNA\tNA\t%d\n", reps_SNP);
	fprintf(stream, "SNP inbreeding coefficient\tNA\tNA\t%d\n", reps_SNP);
        printf("SNP inbreeding coefficient\tNA\tNA\t%d\n", reps_SNP);
}
fprintf(stream, "frequency of major homozygotes\t%6.5f\t%6.5f\tNA\t%d\n", mean_best_freq_major_homo, sd_best_freq_major_homo, niters);
printf("frequency of major homozygotes\t%6.5f\t%6.5f\tNA\t%d\n", mean_best_freq_major_homo, sd_best_freq_major_homo, niters);
fprintf(stream, "frequency of heterozygotes\t%6.5f\t%6.5f\tNA\t%d\n", mean_best_H, sd_best_H, niters);
printf("frequency of heterozygotes\t%6.5f\t%6.5f\tNA\t%d\n", mean_best_H, sd_best_H, niters);
fprintf(stream, "frequency of minor homozygotes\t%6.5f\t%6.5f\tNA\t%d\n", mean_best_freq_minor_homo, sd_best_freq_minor_homo, niters);
printf("frequency of minor homozygotes\t%6.5f\t%6.5f\tNA\t%d\n", mean_best_freq_minor_homo, sd_best_freq_minor_homo, niters);
fprintf(stream, "error rate\t%6.5f\t%6.5f\tNA\t%d\n", mean_best_error, sd_best_error, niters); 
printf("error rate\t%6.5f\t%6.5f\tNA\t%d\n", mean_best_error, sd_best_error, niters);
fprintf(stream, "\n"); 
printf("\n");

// Calculate summary statistics of realized values of parameters in the simulation 
mean_freq_major_homo = (double(sum_freq_major_homo)/double(niters));
mean_freq_minor_homo = (double(sum_freq_minor_homo)/double(niters));
mean_freq_hetero = (double(sum_freq_hetero)/double(niters));
mean_error_rate = (double(sum_error_rate)/double(niters));
mean_major_allele_freq = (double(sum_major_allele_freq)/double(niters));
mean_sD_A = (double(sum_sD_A)/double(niters));
mean_inb_coef = (double(sum_inb_coef)/double(niters));
mean_sq_freq_major_homo = (double(sum_sq_freq_major_homo)/double(niters));
mean_sq_freq_minor_homo = (double(sum_sq_freq_minor_homo)/double(niters));
mean_sq_freq_hetero = (double(sum_sq_freq_hetero)/double(niters));
mean_sq_error_rate = (double(sum_sq_error_rate)/double(niters));
mean_sq_major_allele_freq = (double(sum_sq_major_allele_freq)/double(niters));
mean_sq_sD_A = (double(sum_sq_sD_A)/double(niters));
mean_sq_inb_coef = (double(sum_sq_inb_coef)/double(niters));
var_freq_major_homo = ( niters/(niters-1.0) )*( mean_sq_freq_major_homo - pow(mean_freq_major_homo, 2.0) );
var_freq_minor_homo = ( niters/(niters-1.0) )*( mean_sq_freq_minor_homo - pow(mean_freq_minor_homo, 2.0) );
var_freq_hetero = ( niters/(niters-1.0) )*( mean_sq_freq_hetero - pow(mean_freq_hetero, 2.0) );
var_error_rate = ( niters/(niters-1.0) )*( mean_sq_error_rate - pow(mean_error_rate, 2.0) );
var_major_allele_freq = ( niters/(niters-1.0) )*( mean_sq_major_allele_freq - pow(mean_major_allele_freq, 2.0) );
var_sD_A = ( niters/(niters-1.0) )*( mean_sq_sD_A - pow(mean_sD_A, 2.0) );
var_inb_coef = ( niters/(niters-1.0) )*( mean_sq_inb_coef - pow(mean_inb_coef, 2.0) );
sd_freq_major_homo = pow(var_freq_major_homo, 0.5);
sd_freq_minor_homo = pow(var_freq_minor_homo, 0.5);
sd_freq_hetero = pow(var_freq_hetero, 0.5);
sd_error_rate = pow(var_error_rate, 0.5);
sd_major_allele_freq = pow(var_major_allele_freq, 0.5);
sd_sD_A = pow(var_sD_A, 0.5);
sd_inb_coef = pow(var_inb_coef, 0.5);
if (s_reps_pol > 0) {
	mean_pol_sD_A = sum_pol_sD_A / double(s_reps_pol);
	mean_sq_pol_sD_A = sum_sq_pol_sD_A / double(s_reps_pol);
	var_pol_sD_A = (double(s_reps_pol)/(double(s_reps_pol) - 1.0)) * (mean_sq_pol_sD_A - pow(mean_pol_sD_A,2.0));
	sd_pol_sD_A = pow(var_pol_sD_A,0.5);
	mean_pol_inb_coef = sum_pol_inb_coef / double(s_reps_pol);
	mean_sq_pol_inb_coef = sum_sq_pol_inb_coef / double(s_reps_pol);
	var_pol_inb_coef = (double(s_reps_pol)/(double(s_reps_pol) - 1.0)) * (mean_sq_pol_inb_coef - pow(mean_pol_inb_coef,2.0));
	sd_pol_inb_coef = pow(var_pol_inb_coef,0.5);
}

fprintf(stream, "\n");
printf("\n");
fprintf(stream, "Simulation results\n");
printf("Simulation results\n");
fprintf(stream, "quantity\taverage\tSD\tnumber of data\n");
printf("quantity\taverage\tSD\tnumber of data\n");
fprintf(stream, "frequency of major homozygotes\t%5.4f\t%5.4f\t%d\n", mean_freq_major_homo, sd_freq_major_homo, niters);
printf("frequency of major homozygotes\t%5.4f\t%5.4f\t%d\n", mean_freq_major_homo, sd_freq_major_homo, niters);
fprintf(stream, "frequency of heterozygotes\t%5.4f\t%5.4f\t%d\n", mean_freq_hetero, sd_freq_hetero, niters);
printf( "frequency of heterozygotes\t%5.4f\t%5.4f\t%d\n", mean_freq_hetero, sd_freq_hetero, niters);
fprintf(stream, "frequency of minor homozygotes\t%5.4f\t%5.4f\t%d\n", mean_freq_minor_homo, sd_freq_minor_homo, niters);
printf("frequency of minor homozygotes\t%5.4f\t%5.4f\t%d\n", mean_freq_minor_homo, sd_freq_minor_homo, niters);
fprintf(stream, "error rate\t%5.4f\t%5.4f\t%d\n", mean_error_rate, sd_error_rate, niters);
printf("error rate\t%5.4f\t%5.4f\t%d\n", mean_error_rate, sd_error_rate, niters);
fprintf(stream, "major allele frequency\t%5.4f\t%5.4f\t%d\n", mean_major_allele_freq, sd_major_allele_freq, niters);
printf("major allele frequency\t%5.4f\t%5.4f\t%d\n", mean_major_allele_freq, sd_major_allele_freq, niters);
fprintf(stream, "disequilibrium coefficient\t%5.4f\t%5.4f\t%d\n", mean_sD_A, sd_sD_A, niters);
printf("disequilibrium coefficient\t%5.4f\t%5.4f\t%d\n", mean_sD_A, sd_sD_A, niters);
if (s_reps_pol > 0) {
	fprintf(stream, "polymorphic disequilibrium coefficient\t%6.5f\t%6.5f\t%d\n", mean_pol_sD_A, sd_pol_sD_A, s_reps_pol);
	printf("polymorphic disequilibrium coefficient\t%6.5f\t%6.5f\t%d\n", mean_pol_sD_A, sd_pol_sD_A, s_reps_pol);
} else {
	fprintf(stream, "polymorphic disequilibrium coefficient\tNA\tNA\t%d\n", s_reps_pol);
	printf("polymorphic disequilibrium coefficient\tNA\tNA\t%d\n", s_reps_pol);
}
fprintf(stream, "inbreeding coefficient\t%5.4f\t%5.4f\t%d\n", mean_inb_coef, sd_inb_coef, niters);
printf("inbreeding coefficient\t%5.4f\t%5.4f\t%d\n", mean_inb_coef, sd_inb_coef, niters);
if (s_reps_pol > 0) {
	fprintf(stream, "polymorphic inbreeding coefficient\t%6.5f\t%6.5f\t%d\n", mean_pol_inb_coef, sd_pol_inb_coef, s_reps_pol);
	printf("polymorphic inbreeding coefficient\t%6.5f\t%6.5f\t%d\n", mean_pol_inb_coef, sd_pol_inb_coef, s_reps_pol);
} else {
	fprintf(stream, "polymorphic inbreeding coefficient\tNA\tNA\t%d\n", s_reps_pol);
	printf("polymorphic inbreeding coefficient\tNA\tNA\t%d\n", s_reps_pol);
}

// Close the output file
fclose(stream);

return 0;
}
