// Updated on 09/15/21

/*

Program GFE_v3.0.cpp to estimate allele and genotype frequencies 
from pro files of individual high-throughput sequencing data of 
multiple individuals from a population.  In this version, the
genotype-call (c) mode is added.  An error message is given when
the identity of the minor allele is ambiguous.  The minimum and 
maximum coverage for each individual can be specified by users.
Users can print out the information at all sites in the c mode by 
setting the `-as' option at one.
  

In this version 3.0, the Fixation-indices (F) mode is available.  
The effective number of sampled individuals (Ni) is calculated and 
reported at each site in the F mode.  The error-rate estimate is 
printed out in the F mode.  In addition, the maximum allowed 
error-rate estimate can be specified by users.  The Bayenv (B) mode 
is added to output the allele counts necessary for using the software.  
The site-frequency-spectrum (sfs) mode is available to make bins of 
genotype-frequency estimates equal to the number of individuals with 
read data (Nr). The private-allele mode (p) is added to print out the 
effective number of sampled chromosome (Nc) instead of the effective 
number of sampled individuals (Ni) for the population-structure analysis.  
    
*/

#include <stdlib.h>
#include <iostream>
#include <string>
#include <string.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
using namespace std;

int main(int argc, char *argv[])
{
	// Default values of the options
	const char* in_file_name = {"In_GFE.txt"};
	const char* out_file_name = {"Out_GFE.txt"};
	char mode[10] = {"f"};
	int print_help = 0;
	double cv = 5.991;
	int min_cov = 1;	// minimum coverage required for each individual
	int max_cov = 2000000000;	// maximum coverage allowed for each individual
	int as = 0;	// one if all sites are printed out in the c mode, zero, otherwise.
	double max_e = 1.0;		// maximum allowed error-rate estimate
	int ref_info = 1;
	int sfs = 0;		// one in the sfs mode, zero, otherwise
	const char* pop = {""};

	int argz = 1;	// argument counter

	// Read specified settings
	while( (argz<argc) && (argv[argz][0] == '-') ) {
		if (strcmp(argv[argz], "-h") == 0) {
			print_help = 1;
		} else if (strcmp(argv[argz], "-in") == 0) {
			in_file_name = argv[++argz];
		} else if (strcmp(argv[argz], "-out") == 0) {
			out_file_name = argv[++argz];
		} else if (strcmp(argv[argz], "-mode") == 0) {
			sscanf(argv[++argz], "%s", mode);
		} else if (strcmp(argv[argz], "-cv") == 0) {
			sscanf(argv[++argz], "%lf", &cv);
		} else if (strcmp(argv[argz], "-min_cov") == 0) {
			sscanf(argv[++argz], "%d", &min_cov);
		} else if (strcmp(argv[argz], "-max_cov") == 0) {
			sscanf(argv[++argz], "%d", &max_cov);
		} else if (strcmp(argv[argz], "-as") == 0) {
			sscanf(argv[++argz], "%d", &as);
		} else if (strcmp(argv[argz], "-max_e") == 0) {
			sscanf(argv[++argz], "%lf", &max_e);
		} else if (strcmp(argv[argz], "-pop") == 0) {
			pop = argv[++argz];
		} else if (strcmp(argv[argz], "-ref_info") == 0) {
			sscanf(argv[++argz], "%d", &ref_info);
		} else if (strcmp(argv[argz], "-sfs") == 0) {
			sscanf(argv[++argz], "%d", &sfs);
		} else {
			fprintf(stderr, "unknown option %s\n", argv[argz]);
			print_help = 1;
			break;
		}
		argz++;
	}
	if (print_help) {	// print error/usage message ?
		fprintf(stderr, "USAGE: %s {<options>}\n", argv[0]);
		fprintf(stderr, "	options:\n");
		fprintf(stderr, "	-h: print the usage message\n");
		fprintf(stderr, "	-in <s>: specify the input file name\n");
		fprintf(stderr, "       -out <s>: specify the output file name\n");
		fprintf(stderr, "	-mode <s>: specify the analysis mode\n");
		fprintf(stderr, "       -cv <f>: specify the chi-square critical value for the polymorphism test\n");
		fprintf(stderr, "	-min_cov <d>: specify the minimum coverage required for each individual\n");
		fprintf(stderr, "       -max_cov <d>: specify the maximum coverage allowed for each individual\n");
		fprintf(stderr, "       -as <d>: specify whether all sites are printed out in the c mode\n");
		fprintf(stderr, "       -max_e <f>: specify the maximum allowed error-rate estimate\n");
		fprintf(stderr, "	-pop <s>: specify the population label in the B and F modes\n");
		fprintf(stderr, "       -ref_info <d>: specify whether reference information is shown in the B and F modes\n");
		fprintf(stderr, "       -sfs <d>: specify whether SFS mode is on or not\n");
		exit(1);
	}	
	
	if (cv != 5.991) {
		printf("cv: %f\n", cv);
	}
	if (min_cov != 1) {
		printf("min_cov: %d\n", min_cov);
	}
	if (max_cov != 2000000000) {
		printf("max_cov: %d\n", max_cov);
	}
	if (max_e != 1.0) {
                printf("max_e: %f\n", max_e);
        }

	string line;	// String buffer

	ifstream inputFile(in_file_name);	// Try to open the input file
	if ( !inputFile.is_open() ) {	// Exit on failure
		fprintf(stderr, "Cannot open %s for reading.\n", in_file_name); 
		exit(1);
	}

	// Read the header
	string h_scaf, h_site, h_ref_nuc;	// Stores header names
	getline(inputFile, line);
	istringstream ss(line);
	ss >> h_scaf >> h_site >> h_ref_nuc;
	string str;	// Temporarily stores each individual ID
	vector <string> id_ind;      // Stores individual IDs.
	id_ind.clear();
	while (true) {
		ss >> str;
		id_ind.push_back(str);
		if ( ss.eof() ) {
			break;
		}
	}
	int nsample = id_ind.size();
	printf("%d individuals to be analyzed\n", nsample);

	string s_mode;
	string scaffold;
	int site;
	string ref_nuc;
	string quartet[nsample+1];
	int ag, dg, gg, ig, jg, kg, lg, mg, size_quartet, digit;
	int num_digit[5];
	int count_nuc;
	int read_count[nsample+1][5];
	int ind_coverage[nsample+1];
	double Nc;	// effective number of sampled chromosomes (Mauki and Lynch 2014)
	double Ni;	// effective number of sampled individuals (Mauki and Lynch 2014)
	int Nr;		// Number of individuals with at least one nucleotide read
	int pop_read[5];
	int pop_coverage;
	int n1, n2, n3;
	string major_allele, minor_allele;
	int ml_pop_read[5];
	int ml_ind_read[nsample+1][5];
	int null_ind_read[nsample+1][3];
	double best_error;
	double num_best_p, den_best_p, pre_best_p, best_p, best_q;
	double HWE_maxll;
	double HWE_prob_obs_nuc[nsample+1];
	double null_best_error, null_maxll;
	double null_prob_obs_nuc[nsample+1];
	double size_grid_D_A, mlD_Amin, test_best_p, pre_mlD_Amax, mlD_Amax, mlD_A;
	int nesample;	// effective sample size (equal to Nr in the SFS mode, nsample, otherwise)
	int factor, max_dg;
	double best_error13;
	double prob_geno[4], prob_nuc[4][4];
	double maxll, llhood;
	double prob_obs_nuc[nsample+1];
	double best_D_A, best_f;
	double pre_best_P, pre_best_Q;
	double best_P, best_H, best_Q;
	double best_h;
	double pol_llstat, HWE_llstat;
	int adjust;
	double mlP, mlQ;
	double test_Mac, test_mac;
	int best_Mac, best_mac;	

	FILE *outstream;

	// Open the output file
	outstream = fopen(out_file_name, "w");
	if (outstream == NULL ) {	// Exit on failure
		fprintf(stderr, "Cannot open %s for writing.\n", out_file_name); 
		exit(1);
	}

	s_mode = mode;
	// Print out the field names
	if (s_mode == "f") {
		fprintf(outstream, "scaffold\tsite\tref_nuc\tmajor_allele\tminor_allele\tpop_coverage\tNc\tNr\tbest_Mac\tbest_mac\tbest_p\tbest_q\tbest_error\tnull_best_error\tbest_D_A\tbest_f\tbest_P\tbest_H\tbest_Q\tbest_h\tpol_llstat\tHWE_llstat\n");
		// printf("scaffold\tsite\tref_nuc\tmajor_allele\tminor_allele\tpop_coverage\tNc\tNr\tbest_Mac\tbest_mac\tbest_p\tbest_q\tbest_error\tnull_best_error\tbest_D_A\tbest_f\tbest_P\tbest_H\tbest_Q\tbest_h\tpol_llstat\tHWE_llstat\n");
	} else if (s_mode == "l") {
		fprintf(outstream, "scaffold\tsite\tref_nuc\tn1\tn2\tpop_coverage\tbest_p\tbest_error\tpol_llstat\tHWE_llstat\t");
		// printf("scaffold\tsite\tref_nuc\tn1\tn2\tpop_coverage\tbest_p\tbest_error\tpol_llstat\tHWE_llstat\t");
		for (ig = 0; ig < nsample-1; ig++) {
			fprintf(outstream, "%s\t", id_ind.at(ig).c_str());
			// printf("%s\t", id_ind.at(ig).c_str());
		}
		fprintf(outstream, "%s\n", id_ind.at(ig).c_str());
		// printf("%s\n", id_ind.at(ig).c_str());
	} else if (s_mode == "c") {
		fprintf(outstream, "scaffold\tsite\tref_nuc\tn1\tn2\tpop_coverage\tbest_error\tbest_P\tbest_H\tbest_Q\tpol_llstat\tHWE_llstat\t");
		// printf("scaffold\tsite\tref_nuc\tn1\tn2\tpop_coverage\tbest_error\tbest_P\tbest_H\tbest_Q\tpol_llstat\tHWE_llstat\t");
		for (ig = 0; ig < nsample-1; ig++) {
			fprintf(outstream, "%s\t", id_ind.at(ig).c_str());
			// printf("%s\t", id_ind.at(ig).c_str());
		}
		fprintf(outstream, "%s\n", id_ind.at(ig).c_str());
		// printf("%s\n", id_ind.at(ig).c_str());
	} else if (s_mode == "F") {
		if (ref_info == 1) {
			fprintf(outstream, "scaffold\tsite\tref_nuc\t%s_major_allele\t%s_minor_allele\t%s_pop_coverage\t%s_Ni\t%s_best_Maf\t%s_best_maf\t%s_best_error\t%s_best_H\t%s_pol_llstat\n", pop, pop, pop, pop, pop, pop, pop, pop, pop);
			// printf("scaffold\tsite\tref_nuc\t%s_major_allele\t%s_minor_allele\t%s_pop_coverage\t%s_Ni\t%s_best_Maf\t%s_best_maf\t%s_best_error\t%s_best_H\t%s_pol_llstat\n", pop, pop, pop, pop, pop, pop, pop, pop, pop);
		} else if (ref_info == 0) {
			fprintf(outstream, "%s_major_allele\t%s_minor_allele\t%s_pop_coverage\t%s_Ni\t%s_best_Maf\t%s_best_maf\t%s_best_error\t%s_best_H\t%s_pol_llstat\n", pop, pop, pop, pop, pop, pop, pop, pop, pop);
			// printf("scaffold\tsite\tref_nuc\t%s_major_allele\t%s_minor_allele\t%s_pop_coverage\t%s_Ni\t%s_best_Maf\t%s_best_maf\t%s_best_error, %s_best_H\t%s_pol_llstat\n", pop, pop, pop, pop, pop, pop, pop, pop, pop);
		}
	} else if (s_mode == "B") {
		if (ref_info == 1) {
			fprintf(outstream, "scaffold\tsite\tref_nuc\t%s_major_allele\t%s_minor_allele\t%s_pop_coverage\t%s_Ni\t%s_best_Mac\t%s_best_mac\t%s_best_error\t%s_pol_llstat\n", pop, pop, pop, pop, pop, pop, pop, pop);
			// printf("scaffold\tsite\tref_nuc\t%s_major_allele\t%s_minor_allele\t%s_pop_coverage\t%s_Ni\t%s_best_Mac\t%s_best_mac\t%s_best_error\t%s_pol_llstat\n", pop, pop, pop, pop, pop, pop, pop, pop);
		} else if (ref_info == 0) {
			fprintf(outstream, "%s_major_allele\t%s_minor_allele\t%s_pop_coverage\t%s_Ni\t%s_best_Mac\t%s_best_mac\t%s_best_error\t%s_pol_llstat\n", pop, pop, pop, pop, pop, pop, pop, pop);
			// printf("%s_major_allele\t%s_minor_allele\t%s_pop_coverage\t%s_Ni\t%s_best_Mac\t%s_best_mac\t%s_best_error\t%s_pol_llstat\n", pop, pop, pop, pop, pop, pop, pop, pop);
		}
	} else if (s_mode == "p") {
                if (ref_info == 1) {
                        fprintf(outstream, "scaffold\tsite\tref_nuc\t%s_major_allele\t%s_minor_allele\t%s_pop_coverage\t%s_Nc\t%s_best_Maf\t%s_best_maf\t%s_best_error\t%s_best_H\t%s_pol_llstat\n", pop, pop, pop, pop, pop, pop, pop, pop, pop);
                        // printf("scaffold\tsite\tref_nuc\t%s_major_allele\t%s_minor_allele\t%s_pop_coverage\t%s_Nc\t%s_best_Maf\t%s_best_maf\t%s_best_error\t%s_best_H\t%s_pol_llstat\n", pop, pop, pop, pop, pop, pop, pop, pop, pop);
                } else if (ref_info == 0) {
                        fprintf(outstream, "%s_major_allele\t%s_minor_allele\t%s_pop_coverage\t%s_Nc\t%s_best_Maf\t%s_best_maf\t%s_best_error\t%s_best_H\t%s_pol_llstat\n", pop, pop, pop, pop, pop, pop, pop, pop, pop);
                        // printf("scaffold\tsite\tref_nuc\t%s_major_allele\t%s_minor_allele\t%s_pop_coverage\t%s_Nc\t%s_best_Maf\t%s_best_maf\t%s_best_error, %s_best_H\t%s_pol_llstat\n", pop, pop, pop, pop, pop, pop, pop, pop, pop);
                }
	}
	// Read the main data	
	while ( getline(inputFile, line) ) {
		istringstream ss(line);
		ss >> scaffold >> site >> ref_nuc;
		for (ig = 1; ig <= nsample; ig++) {
			ss >> quartet[ig];
		}
		for (jg = 1; jg <= 4; jg++) {
                        pop_read[jg] = 0;
                }
		pop_coverage = 0;
		Nc = 0.0;
		Ni = 0.0;
		Nr = 0; 
		for (ig = 1; ig <= nsample; ig++) {
			size_quartet = quartet[ig].size();
			jg = 1;
			digit = 0;
			ind_coverage[ig] = 0;
			for (kg = 0; kg < size_quartet; kg++) {
				if ( quartet[ig].at(kg) == '/') {
					mg = pow(10,digit-1);
					count_nuc = 0;
					lg = 0;
					while (lg <= digit-1) {
						count_nuc = count_nuc + mg*num_digit[lg];
						mg = mg/10;
						num_digit[lg] = 0;
						lg = lg + 1;
					}
					read_count[ig][jg] = count_nuc;
					pop_read[jg] = pop_read[jg] + read_count[ig][jg];
					ind_coverage[ig] = ind_coverage[ig] + count_nuc;
					jg = jg + 1;
					count_nuc = 0;
					digit = 0;
				} else {
					num_digit[digit] = quartet[ig].at(kg) - '0';
					digit = digit + 1;
				}
				if (kg == size_quartet - 1) {
					mg = pow(10,digit-1);
					count_nuc = 0;
					lg = 0;
					while (lg <= digit-1) {
						count_nuc = count_nuc + mg*num_digit[lg];
						mg = mg/10;
						num_digit[lg] = 0;
						lg = lg + 1;
					}
					read_count[ig][jg] = count_nuc;
					ind_coverage[ig] = ind_coverage[ig] + count_nuc;
					pop_read[jg] = pop_read[jg] + read_count[ig][jg];
					jg = 1;
                		}			 
			}
			if (ind_coverage[ig] < min_cov || ind_coverage[ig] > max_cov) {
				for (jg = 1; jg <= 4; jg++) {
					read_count[ig][jg] = 0;
				}
				ind_coverage[ig] = 0;
			}
			if (ind_coverage[ig] > 0) {
				Nc = Nc + (double)2.0 - pow((double)0.5,ind_coverage[ig]-1);
				Ni = Ni + (double)1.0 - pow((double)0.5,ind_coverage[ig]);
				Nr = Nr + 1;
			}	
		}
		for (jg = 1; jg <= 4; jg++) {
                	pop_coverage = pop_coverage + pop_read[jg];
                }
		if (sfs == 1) {
			nesample = Nr;
		} else {
			nesample = nsample;
		}
		// ML estimation starts here
		if (pop_coverage > 0) {
			// Find the two most abundant nucleotides at the site and consider them as candidate alleles at the site.
			if (pop_read[1] >= pop_read[2]) {
				n1 = 1;
				n2 = 2;
			} else {
				n1 = 2;
				n2 = 1;
			}
			if (pop_read[3] > pop_read[n1]) {
				n3 = n2;
				n2 = n1;
				n1 = 3;
			} else if (pop_read[3] > pop_read[n2]) {
				n3 = n2;
				n2 = 3;
			} else {
				n3 = 3;
			}	
			if (pop_read[4] > pop_read[n1]) {
				n3 = n2;
				n2 = n1;
				n1 = 4;
			} else if (pop_read[4] > pop_read[n2]) {
				n3 = n2;
				n2 = 4;
			} else if (pop_read[4] > pop_read[n3]) {
				n3 = 4;
			}
			if (n1 == 1) {
				major_allele = "A";
			} else if (n1 == 2) {
				major_allele = "C";
			} else if (n1 == 3) {
				major_allele = "G";
			} else if (n1 == 4) {
				major_allele = "T";
			}
			if (n2 == 1) {
                        	minor_allele = "A";
                	} else if (n2 == 2) {
                        	minor_allele = "C";
                	} else if (n2 == 3) {
                        	minor_allele = "G";
                	} else if (n2 == 4) {
                        	minor_allele = "T";
                	}
			// Carry out the ML estimation only when there are at least two different nucleotides in the population sample
			if (pop_read[n1] < pop_coverage) {
				// Count the number of different types of reads at the site for use in the ML estimation
				ml_pop_read[1] = pop_read[n1];
				ml_pop_read[2] = pop_read[n2];
				ml_pop_read[3] = pop_coverage - ml_pop_read[1] - ml_pop_read[2];
				for (ig = 1; ig <= nsample; ig++) {
					ml_ind_read[ig][1] = read_count[ig][n1];
					ml_ind_read[ig][2] = read_count[ig][n2];
					ml_ind_read[ig][3] = ind_coverage[ig] - ml_ind_read[ig][1] - ml_ind_read[ig][2];
					null_ind_read[ig][1] = read_count[ig][n1];
					null_ind_read[ig][2] = ind_coverage[ig] - null_ind_read[ig][1];
				}
				num_best_p = 2.0*ml_pop_read[1] - ml_pop_read[3];
				den_best_p = 2.0*(ml_pop_read[1] + ml_pop_read[2] - ml_pop_read[3]);
				if (den_best_p == 0.0) {
					fprintf(stderr, "ML estimates not found at site %d on %s\n", site, scaffold.c_str());
					// Print out the results in the f, B, and F modes
					if (s_mode == "f") {
						fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
						// printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
					} else if (s_mode == "c" && as == 1) {
						fprintf(outstream, "%s\t%d\t%s\t%d\t%d\t%d\tNA\tNA\tNA\tNA\tNA\tNA\t", scaffold.c_str(), site, ref_nuc.c_str(), n1, n2, pop_coverage);
						// printf("%s\t%d\t%s\t%d\t%d\t%d\tNA\tNA\tNA\tNA\tNA\tNA\t", scaffold.c_str(), site, ref_nuc.c_str(), n1, n2, pop_coverage);
						for (ig=1;ig<nsample;ig++) {
							fprintf(outstream, "%s\t", quartet[ig].c_str());
							// printf("%s\t", quartet[ig].c_str());
						}
						fprintf(outstream, "%s\n", quartet[nsample].c_str());
						// printf("%s\n", quartet[nsample].c_str());
					} else if (s_mode == "F") {
						if (ref_info == 1) {
							fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
							// printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
						} else if (ref_info == 0) {
							fprintf(outstream, "NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
							// printf("NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
						}
					} else if (s_mode == "B") {
						if (ref_info == 1) {
							fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
							// printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
						} else if (ref_info == 0) {
							fprintf(outstream, "NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
							// printf("NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
						}
					} else if (s_mode == "p") {
                                                if (ref_info == 1) {
                                                        fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
                                                        // printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
                                                } else if (ref_info == 0) {
                                                        fprintf(outstream, "NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
                                                        // printf("NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
                                                }
                                        }
				} else {
					best_error = 1.5*(ml_pop_read[3]/(double)pop_coverage);
					pre_best_p = num_best_p/den_best_p;
					if ( 2.0*nesample*pre_best_p - (int)(2.0*nesample*pre_best_p) >= 0.5 ) {
						best_p = (double)(1.0/(2.0*nesample))*( (int)(2.0*nesample*pre_best_p)+1 );
					} else {
						best_p = (double)(1.0/(2.0*nesample))*( (int)(2.0*nesample*pre_best_p) );
					}
					if (best_p > 1.0 - 1.0e-8) {
						best_p = (double)(1.0/(2.0*nesample))*( (int)(2.0*nesample*pre_best_p) );
					}
					// Calculate the corresponding maximum log-likelihood under the null hypothesis of monomorphism.
					null_best_error = (double)(pop_coverage - ml_pop_read[1])/(double)pop_coverage;
					// Sum the log-likelihoods under the null hypothesis of monomorphism
					null_maxll = 0.0;
					for (ig = 1; ig <= nsample; ig++) {
						if (ind_coverage[ig] > 0) {
							null_prob_obs_nuc[ig] = pow( 1.0-null_best_error, (double)null_ind_read[ig][1] )*pow(null_best_error/3.0, (double)(double)null_ind_read[ig][2]);
							if (null_prob_obs_nuc[ig] > 0.0) {
								null_maxll = null_maxll + log(null_prob_obs_nuc[ig]);
							} else {
								null_maxll = -10000000001.0;
							}
						}
					}
					// Estimate the disequilibrium coefficient by a grid search
					size_grid_D_A = 1.0/(double)nesample;
					if ( -pow(best_p,2.0) >= -pow(1.0-best_p,2.0) ) {
						mlD_Amin = -pow(best_p,2.0);
					} else {
						mlD_Amin = -pow(1.0-best_p,2.0);
					}
					test_best_p = best_p*(double)nesample-(int)(best_p*nesample);
					if (test_best_p < 1.0e-8) {
						mlD_Amax = best_p*(1.0-best_p);
					} else {
						pre_mlD_Amax = best_p*(1.0-best_p);
						factor = (int)( (pre_mlD_Amax-mlD_Amin)/size_grid_D_A );
						mlD_Amax = mlD_Amin + size_grid_D_A*factor;
					}
					mlD_A = mlD_Amin - size_grid_D_A;
					max_dg = (int)( nesample*(mlD_Amax-mlD_Amin) ) + 2;
					maxll = -10000000000.0;
					// Loop over the candidate disequilibrium coefficients at the site
					for (dg = 1; dg <= max_dg; ++dg) {
						if (dg == max_dg) {
							mlD_A = mlD_Amax;
						} else {
							mlD_A = mlD_A + size_grid_D_A;
						}
						// Calculate the probabilities for use in the likelihood function
						best_error13 = best_error/3.0;
						// Calculate the probability of each of the three genotypes for an individual
						prob_geno[1] = pow(best_p,2.0) + mlD_A;	 
						prob_geno[2] = 2.0*( best_p*(1.0-best_p)-mlD_A );
						prob_geno[3] = pow(1.0-best_p,2.0) + mlD_A;
						// Calculate the probability of an observed nucleotide given a genotype of the individual at the site
						prob_nuc[1][1] = 1.0 - best_error;
						prob_nuc[1][2] = best_error13;
						prob_nuc[1][3] = best_error13;
						prob_nuc[2][1] = 0.5 - best_error13;
						prob_nuc[2][2] = 0.5 - best_error13;
						prob_nuc[2][3] = best_error13;
						prob_nuc[3][1] = best_error13;
						prob_nuc[3][2] = 1.0 - best_error;
						prob_nuc[3][3] = best_error13;
						// Sum the log-likelihoods over the individuals
						llhood = 0.0;
						for (ig = 1; ig <= nsample; ig++) {
							if (ind_coverage[ig] > 0) {	
								// Sum the probabilities over the genotypes of the individual
								prob_obs_nuc[ig] = 0.0;
								for (gg = 1; gg <= 3; ++gg) {
									prob_obs_nuc[ig] = prob_obs_nuc[ig] + prob_geno[gg]*pow(prob_nuc[gg][1], (double)ml_ind_read[ig][1])*pow(prob_nuc[gg][2], (double)ml_ind_read[ig][2])*pow(prob_nuc[gg][3], (double)ml_ind_read[ig][3]);
								}
								if (prob_obs_nuc[ig] > 0.0) {
									llhood = llhood + log(prob_obs_nuc[ig]);
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
					} // End the loop for disequilibrium coefficients at the site
					pre_best_P = pow(best_p,2.0) + best_D_A;
					pre_best_Q = pow(1.0-best_p,2.0) + best_D_A;
					if (pre_best_P < 0.0) {
						best_P = 0.0;
					} else {
						if ( nesample*pre_best_P - (int)(nesample*pre_best_P) >= 0.5 ) {
							best_P = (double)(1.0/(double)nesample)*( (int)(nesample*pre_best_P)+1 );
						} else {
							best_P = (double)(1.0/(double)nesample)*( (int)(nesample*pre_best_P) );
						}
					}
					if (pre_best_Q < 0.0) {
						best_Q = 0.0;
					} else {
						if ( nesample*pre_best_Q - (int)(nesample*pre_best_Q) >= 0.5 ) {
							best_Q = (double)(1.0/(double)nesample)*( (int)(nesample*pre_best_Q)+1 );
						} else {
							best_Q = (double)(1.0/(double)nesample)*( (int)(nesample*pre_best_Q) );
						}
					}

					// The adjustment of the ML genotype-frequency estimates starts here
					// Adjust the ML estimates by a grid search
					adjust = 1;
					while (adjust == 1) {
						adjust = 0;
						for (ag=1; ag<=8; ++ag) {
							if (ag == 1) {
								mlP = best_P - size_grid_D_A;
								mlQ = best_Q - size_grid_D_A;
							} else if (ag == 2) {
								mlP = best_P - size_grid_D_A;
								mlQ = best_Q;
							} else if (ag == 3) {
								mlP = best_P - size_grid_D_A;
								mlQ = best_Q + size_grid_D_A;
							} else if (ag == 4) {
								mlP = best_P;
								mlQ = best_Q - size_grid_D_A;
							} else if (ag == 5) {
								mlP = best_P;
								mlQ = best_Q + size_grid_D_A;
							} else if (ag == 6) {
								mlP = best_P + size_grid_D_A;
								mlQ = best_Q - size_grid_D_A;
							} else if (ag == 7) {
								mlP = best_P + size_grid_D_A;
								mlQ = best_Q;
							} else if (ag == 8) {
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
								for (ig = 1; ig <= nsample; ig++) {
									if (ind_coverage[ig] > 0) {
										// Sum the probabilities over the genotypes of the individual
										prob_obs_nuc[ig] = 0.0;
										for (gg = 1; gg <= 3; ++gg) {
											prob_obs_nuc[ig] = prob_obs_nuc[ig] + prob_geno[gg]*pow(prob_nuc[gg][1], (double)ml_ind_read[ig][1])*pow(prob_nuc[gg][2], (double)ml_ind_read[ig][2])*pow(prob_nuc[gg][3], (double)ml_ind_read[ig][3]);
										}
										if (prob_obs_nuc[ig] > 0.0) {
											llhood = llhood + log(prob_obs_nuc[ig]);
										} else {
											llhood = -10000000001.0;
										}
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
					if (maxll > -10000000000.0) {
						// Calculate the maximum log-likelihood under the hypothesis of HWE
						best_p = best_P + ( (double)1.0-best_P-best_Q )/(double)2.0;
						HWE_maxll = 0.0;
						for (ig = 1; ig <= nsample; ig++) {
							if (ind_coverage[ig] > 0) {
								HWE_prob_obs_nuc[ig] = pow( best_p, 2.0 )*pow( 1.0-best_error, (double)ml_ind_read[ig][1] )*pow( best_error/3.0, (double)ml_ind_read[ig][2] )*pow( best_error/3.0, (double)ml_ind_read[ig][3] ) + 2.0*best_p*(1.0-best_p)*pow( 0.5-best_error/3.0, (double)(ml_ind_read[ig][1]+ml_ind_read[ig][2]) )*pow( best_error/3.0, (double)ml_ind_read[ig][3] ) + pow( 1.0-best_p, 2.0 )*pow( best_error/3.0, (double)ml_ind_read[ig][1] )*pow( 1.0-best_error, (double)ml_ind_read[ig][2] )*pow( best_error/3.0, (double)ml_ind_read[ig][3] );
								if (HWE_prob_obs_nuc[ig] > 0.0) {
									HWE_maxll = HWE_maxll + log(HWE_prob_obs_nuc[ig]);
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
						// The adjustment of the ML estimates ends here
						best_H = (double)1.0 - best_P - best_Q;
						// Fix tiny numerical errors
						if (best_H < 0.0) {
							best_H = 0.0;
						}
						best_p = best_P + best_H/(double)2.0;
						best_q = 1.0 - best_p;
						best_D_A = best_P - pow(best_p,2.0);
						best_h = 2.0*best_p*(1.0-best_p);
						pol_llstat = 2.0*(maxll - null_maxll);
						HWE_llstat = 2.0*(maxll - HWE_maxll);
						if (best_p < 1.0) {
							if (pop_read[n2] == pop_read[n3]) {
                                				// printf("pop_read[1]: %d\tpop_read[2]: %d\tpop_read[3]: %d\tpop_read[4]: %d\n", pop_read[1], pop_read[2], pop_read[3], pop_read[4]);
								// printf("n1: %d\tn2: %d\tn3: %d\n", n1, n2, n3);
                                				// printf("pop_read[n1]: %d\tpop_read[n2]: %d\tpop_read[n3]: %d\n", pop_read[n1], pop_read[n2], pop_read[n3]);
                                				fprintf(stderr, "The identity of the minor allele is ambiguous at site %d on %s\n", site, scaffold.c_str());
                        				}
							best_f = (best_h - best_H)/best_h;
							if (best_p >= best_q) {
								if (s_mode == "f") {
									test_Mac = 2.0*nesample*best_p - (int)(2.0*nesample*best_p);
                                                                	if ( test_Mac >= 0.5) {
                                                                        	best_Mac = (int)(2.0*nesample*best_p) + 1;
                                                                	} else {
                                                                        	best_Mac = (int)(2.0*nesample*best_p);
                                                                	}
                                                                	test_mac = 2.0*nesample*best_q - (int)(2.0*nesample*best_q);
                                                                	if ( test_mac > 0.5) {
                                                                        	best_mac = (int)(2.0*nesample*best_q) + 1;
                                                                	} else {
                                                                        	best_mac = (int)(2.0*nesample*best_q);
                                                                	}
									if (best_error <= max_e) {
										fprintf(outstream, "%s\t%d\t%s\t%s\t%s\t%d\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), minor_allele.c_str(), pop_coverage, Nc, Nr, best_Mac, best_mac, best_p, best_q, best_error, null_best_error, best_D_A, best_f, best_P, best_H, best_Q, best_h, pol_llstat, HWE_llstat);
										// printf("%s\t%d\t%s\t%s\t%s\t%d\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), minor_allele.c_str(), pop_coverage, Nc, Nr, best_Mac, best_mac, best_p, best_q, best_error, null_best_error, best_D_A, best_f, best_P, best_H, best_Q, best_h, pol_llstat, HWE_llstat);
									} else {
										fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
										// printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
									}
								} else if (s_mode == "l") {
									// Print out the results only when the site is significantly polymorphic at the 5% level
									if ( pol_llstat > cv && best_error <= max_e) {
										fprintf(outstream, "%s\t%d\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), n1, n2, pop_coverage, best_p, best_error, pol_llstat, HWE_llstat);
										// printf("%s\t%d\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), n1, n2, pop_coverage, best_p, best_error, pol_llstat, HWE_llstat);
										for (ig=1;ig<nsample;ig++) {
											fprintf(outstream, "%s\t", quartet[ig].c_str());
											// printf("%s\t", quartet[ig].c_str());
										}
										fprintf(outstream, "%s\n", quartet[nsample].c_str());
										// printf("%s\n", quartet[nsample].c_str());
									}
								} else if (s_mode == "c") {
									if (as == 1) {
										fprintf(outstream, "%s\t%d\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), n1, n2, pop_coverage, best_error, best_P, best_H, best_Q, pol_llstat, HWE_llstat);
										// printf("%s\t%d\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), n1, n2, pop_coverage, best_error, best_P, best_H, best_Q, pol_llstat, HWE_llstat);
										for (ig=1;ig<nsample;ig++) {
											fprintf(outstream, "%s\t", quartet[ig].c_str());
											// printf("%s\t", quartet[ig].c_str());
										}
										fprintf(outstream, "%s\n", quartet[nsample].c_str());
										// printf("%s\n", quartet[nsample].c_str());
									} else {
										// Print out the results only when the site is significantly polymorphic at the 5% level
										if ( pol_llstat > cv && best_error <= max_e) {
											fprintf(outstream, "%s\t%d\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), n1, n2, pop_coverage, best_error, best_P, best_H, best_Q, pol_llstat, HWE_llstat);
											// printf("%s\t%d\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), n1, n2, pop_coverage, best_error, best_P, best_H, best_Q, pol_llstat, HWE_llstat);
											for (ig=1;ig<nsample;ig++) {
												fprintf(outstream, "%s\t", quartet[ig].c_str());
												// printf("%s\t", quartet[ig].c_str());
											}
											fprintf(outstream, "%s\n", quartet[nsample].c_str());
											// printf("%s\n", quartet[nsample].c_str());
										}
									}
								} else if (s_mode == "F") {
									if (ref_info == 1) {
										if (best_error <= max_e) {
											fprintf(outstream, "%s\t%d\t%s\t%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), minor_allele.c_str(), pop_coverage, Ni, best_p, best_q, best_error, best_H, pol_llstat);
											// printf("%s\t%d\t%s\t%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), minor_allele.c_str(), pop_coverage, Ni, best_p, best_q, best_error, best_H, pol_llstat);
										} else {
											fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
											// printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
										}
									} else if (ref_info == 0) {
										if (best_error <= max_e) {
											fprintf(outstream, "%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", major_allele.c_str(), minor_allele.c_str(), pop_coverage, Ni, best_p, best_q, best_error, best_H, pol_llstat);
											// printf("%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", major_allele.c_str(), minor_allele.c_str(), pop_coverage, Ni, best_p, best_q, best_error, best_H, pol_llstat);
										} else {
											fprintf(outstream, "NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
											// printf("NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
										}
									} 
								} else if (s_mode == "B") {
									test_Mac = 2.0*nesample*best_p - (int)(2.0*nesample*best_p);
                                                                        if ( test_Mac >= 0.5) {
                                                                                best_Mac = (int)(2.0*nesample*best_p) + 1;
                                                                        } else {
                                                                                best_Mac = (int)(2.0*nesample*best_p);
                                                                        }
                                                                        test_mac = 2.0*nesample*best_q - (int)(2.0*nesample*best_q);
                                                                        if ( test_mac > 0.5) {
                                                                                best_mac = (int)(2.0*nesample*best_q) + 1;
                                                                        } else {
                                                                                best_mac = (int)(2.0*nesample*best_q);
                                                                        }
									if (ref_info == 1) {
										if (best_error <= max_e) {
											fprintf(outstream, "%s\t%d\t%s\t%s\t%s\t%d\t%f\t%d\t%d\t%f\t%f\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), minor_allele.c_str(), pop_coverage, Ni, best_Mac, best_mac, best_error, pol_llstat);
											// printf("%s\t%d\t%s\t%s\t%s\t%d\t%f\t%d\t%d\t%f\t%f\t%f\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), minor_allele.c_str(), pop_coverage, Ni, best_Mac, best_mac, best_error, pol_llstat);
										} else {
											fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
											// printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\n");
										}
									} else if (ref_info == 0) {
										if (best_error <= max_e) {
											fprintf(outstream, "%s\t%s\t%d\t%f\t%d\t%d\t%f\t%f\n", major_allele.c_str(), minor_allele.c_str(), pop_coverage, Ni, best_Mac, best_mac, best_error, pol_llstat);
											// printf("%s\t%s\t%d\t%f\t%d\t%d\t%f\t%f\n", major_allele.c_str(), minor_allele.c_str(), pop_coverage, Ni, best_Mac, best_mac, best_error, pol_llstat);
										} else {
											fprintf(outstream, "NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
											// printf("NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
										}
									}
								} else if (s_mode == "p") {
                                                                        if (ref_info == 1) {
                                                                                if (best_error <= max_e) {
                                                                                        fprintf(outstream, "%s\t%d\t%s\t%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), minor_allele.c_str(), pop_coverage, Nc, best_p, best_q, best_error, best_H, pol_llstat);
                                                                                        // printf("%s\t%d\t%s\t%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), minor_allele.c_str(), pop_coverage, Nc, best_p, best_q, best_error, best_H, pol_llstat);
                                                                                } else {
                                                                                        fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
                                                                                        // printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
                                                                                }
                                                                        } else if (ref_info == 0) {
                                                                                if (best_error <= max_e) {
                                                                                        fprintf(outstream, "%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", major_allele.c_str(), minor_allele.c_str(), pop_coverage, Nc, best_p, best_q, best_error, best_H, pol_llstat);
                                                                                        // printf("%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", major_allele.c_str(), minor_allele.c_str(), pop_coverage, Nc, best_p, best_q, best_error, best_H, pol_llstat);
                                                                                } else {
                                                                                        fprintf(outstream, "NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
                                                                                        // printf("NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
                                                                                }
                                                                        }
                                                                }   
							} else {
								if (s_mode == "f") {
									test_Mac = 2.0*nesample*best_q - (int)(2.0*nesample*best_q);
									if ( test_Mac >= 0.5) {
                                                                		best_Mac = (int)(2.0*nesample*best_q) + 1;
                                                        		} else {
                                                                		best_Mac = (int)(2.0*nesample*best_q);
                                                        		}
									test_mac = 2.0*nesample*best_p - (int)(2.0*nesample*best_p);
                                                        		if ( test_mac >= 0.5) {
                                                                		best_mac = (int)(2.0*nesample*best_p) + 1;
                                                        		} else {
                                                                		best_mac = (int)2.0*(nesample*best_p);
                                                        		}
									if (best_error <= max_e) {
										fprintf(outstream, "%s\t%d\t%s\t%s\t%s\t%d\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", scaffold.c_str(), site, ref_nuc.c_str(), minor_allele.c_str(), major_allele.c_str(), pop_coverage, Nc, Nr, best_Mac, best_mac, best_q, best_p, best_error, null_best_error, best_D_A, best_f, best_Q, best_H, best_P, best_h, pol_llstat, HWE_llstat);
										// printf("%s\t%d\t%s\t%s\t%s\t%d\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", scaffold.c_str(), site, ref_nuc.c_str(), minor_allele.c_str(), major_allele.c_str(), pop_coverage, Nc, Nr, best_Mac, best_mac, best_q, best_p, best_error, null_best_error, best_D_A, best_f, best_Q, best_H, best_P, best_h, pol_llstat, HWE_llstat);
									} else {
										fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
										// printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
									}
								} else if (s_mode == "l") {
									// Print out the results only when the site is significantly polymorphic at the 5% level
									if ( pol_llstat > cv && best_error <= max_e) {
										fprintf(outstream, "%s\t%d\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), n2, n1, pop_coverage, best_q, best_error, pol_llstat, HWE_llstat);
										// printf("%s\t%d\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), n2, n1, pop_coverage, best_q, best_error, pol_llstat, HWE_llstat);
										for (ig=1;ig<nsample;ig++){
											fprintf(outstream, "%s\t", quartet[ig].c_str());
											// printf("%s\t", quartet[ig].c_str());
										}
										fprintf(outstream, "%s\n", quartet[nsample].c_str());
										// printf("%s\n", quartet[nsample].c_str());
									}
								} else if (s_mode == "c") {
									if (as == 1) {
										fprintf(outstream, "%s\t%d\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), n2, n1, pop_coverage, best_error, best_Q, best_H, best_P, pol_llstat, HWE_llstat);
										// printf("%s\t%d\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), n2, n1, pop_coverage, best_error, best_Q, best_H, best_P, pol_llstat, HWE_llstat);
										for (ig=1;ig<nsample;ig++) {
											fprintf(outstream, "%s\t", quartet[ig].c_str());
											// printf("%s\t", quartet[ig].c_str());
										}
										fprintf(outstream, "%s\n", quartet[nsample].c_str());
										// printf("%s\n", quartet[nsample].c_str());
									} else {
                                                                        	// Print out the results only when the site is significantly polymorphic at the 5% level
                                                                        	if ( pol_llstat > cv && best_error <= max_e) {
 											fprintf(outstream, "%s\t%d\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), n2, n1, pop_coverage, best_error, best_Q, best_H, best_P, pol_llstat, HWE_llstat);
											// printf("%s\t%d\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), n2, n1, pop_coverage, best_error, best_Q, best_H, best_P, pol_llstat, HWE_llstat);
											for (ig=1;ig<nsample;ig++) {
												fprintf(outstream, "%s\t", quartet[ig].c_str());
												// printf("%s\t", quartet[ig].c_str());
											}
											fprintf(outstream, "%s\n", quartet[nsample].c_str());
											// printf("%s\n", quartet[SAMPLE_N].c_str());
										}
									}
								} else if (s_mode == "F") {
									if (ref_info == 1) {
										if (best_error <= max_e) {
											fprintf(outstream, "%s\t%d\t%s\t%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", scaffold.c_str(), site, ref_nuc.c_str(), minor_allele.c_str(), major_allele.c_str(), pop_coverage, Ni, best_q, best_p, best_error, best_H, pol_llstat);
											// printf("%s\t%d\t%s\t%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", scaffold.c_str(), site, ref_nuc.c_str(), minor_allele.c_str(), major_allele.c_str(), pop_coverage, Ni, best_q, best_p, best_error, best_H, pol_llstat);
										} else {
											fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
											// printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
										}
									} else if (ref_info == 0) {
										if (best_error <= max_e) {
											fprintf(outstream, "%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", minor_allele.c_str(), major_allele.c_str(), pop_coverage, Ni, best_q, best_p, best_error, best_H, pol_llstat);
											// printf("%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", minor_allele.c_str(), major_allele.c_str(), pop_coverage, Ni, best_q, best_p, best_error, best_H, pol_llstat);
										} else {
											fprintf(outstream, "NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
											// printf("NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
										}
									} 
								} else if (s_mode == "B") {
									test_Mac = 2.0*nesample*best_q - (int)(2.0*nesample*best_q);
                                                                        if ( test_Mac >= 0.5) {
                                                                                best_Mac = (int)(2.0*nesample*best_q) + 1;
                                                                        } else {
                                                                                best_Mac = (int)(2.0*nesample*best_q);
                                                                        }
                                                                        test_mac = 2.0*nesample*best_p - (int)(2.0*nesample*best_p);
                                                                        if ( test_mac >= 0.5) {
                                                                                best_mac = (int)(2.0*nesample*best_p) + 1;
                                                                        } else {
                                                                                best_mac = (int)2.0*(nesample*best_p);
                                                                        }
									if (ref_info == 1) {
										if (best_error <= max_e) {
											fprintf(outstream, "%s\t%d\t%s\t%s\t%s\t%d\t%f\t%d\t%d\t%f\t%f\n", scaffold.c_str(), site, ref_nuc.c_str(), minor_allele.c_str(), major_allele.c_str(), pop_coverage, Ni, best_Mac, best_mac, best_error, pol_llstat);
											// printf("%s\t%d\t%s\t%s\t%s\t%d\t%f\t%d\t%d\t%f\t%f\n", scaffold.c_str(), site, ref_nuc.c_str(), minor_allele.c_str(), major_allele.c_str(), $pop_coverage, Ni, best_Mac, best_mac, best_error, pol_llstat); 
										} else {
											fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
											// printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
										}
									} else if (ref_info == 0) {
										if (best_error <= max_e) {
											fprintf(outstream, "%s\t%s\t%d\t%f\t%d\t%d\t%f\t%f\n", minor_allele.c_str(), major_allele.c_str(), pop_coverage, Ni, best_Mac, best_mac, best_error, pol_llstat);
											// printf("%s\t%s\t%d\t%f\t%d\t%d\t%f\t%f\n", minor_allele.c_str(), major_allele.c_str(), pop_coverage, Ni, best_Mac, best_mac, best_error, pol_llstat);
										} else {
											fprintf(outstream, "NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
											// printf("NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
										}
									}
								} else if (s_mode == "p") {
									if (ref_info == 1) {
										if (best_error <= max_e) {
											fprintf(outstream, "%s\t%d\t%s\t%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", scaffold.c_str(), site, ref_nuc.c_str(), minor_allele.c_str(), major_allele.c_str(), pop_coverage, Nc, best_q, best_p, best_error, best_H, pol_llstat);
											// printf("%s\t%d\t%s\t%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", scaffold.c_str(), site, ref_nuc.c_str(), minor_allele.c_str(), major_allele.c_str(), pop_coverage, Nc, best_q, best_p, best_error, best_H, pol_llstat);
										} else {
											fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
											// printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
										}
									} else if (ref_info == 0) {
										if (best_error <= max_e) {
											fprintf(outstream, "%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", minor_allele.c_str(), major_allele.c_str(), pop_coverage, Nc, best_q, best_p, best_error, best_H, pol_llstat);
											// printf("%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n", minor_allele.c_str(), major_allele.c_str(), pop_coverage, Nc, best_q, best_p, best_error, best_H, pol_llstat);
										} else {
											fprintf(outstream, "NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
											// printf("NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
										}
									}
								} 
							}				
						} else {
							// Print out the monomorphism results for the f, c (with the as option equal to one), B, F, and p modes
							if (s_mode == "f") { 
								test_Mac = 2.0*nesample*best_p - (int)(2.0*nesample*best_p);
								if ( test_Mac >= 0.5) {
                                                			best_Mac = (int)(2.0*nesample*best_p) + 1;
                                                		} else {
                                                        		best_Mac = (int)(2.0*nesample*best_p);
                                                		}
								test_mac = 2.0*nesample*best_q - (int)(2.0*nesample*best_q);
                                                		if ( test_mac >= 0.5) {
                                                        		best_mac = (int)(2.0*nesample*best_q) + 1;
                                                		} else {
                                                        		best_mac = (int)(2.0*nesample*best_q);
                                                		}
								if (best_error <= max_e) {
									fprintf(outstream, "%s\t%d\t%s\t%s\tNA\t%d\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\tNA\tNA\t%f\t%f\t%f\t%f\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), pop_coverage, Nc, Nr, best_Mac, best_mac, best_p, best_q, best_error, null_best_error, best_P, best_H, best_Q, best_h);
									// printf("%s\t%d\t%s\t%s\tNA\t%d\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\tNA\tNA\t%f\t%f\t%f\t%f\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), pop_coverage, Nc, Nr, best_Mac, best_mac, best_p, best_q, best_error, null_best_error, best_P, best_H, best_Q, best_h);
								} else {
									fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
									// printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
								}
							} else if (s_mode == "c" && as == 1) {
								fprintf(outstream, "%s\t%d\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), n1, n2, pop_coverage, best_error, best_P, best_H, best_Q, pol_llstat, HWE_llstat);
								// printf("%s\t%d\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t", scaffold.c_str(), site, ref_nuc.c_str(), n1, n2, pop_coverage, best_error, best_P, best_H, best_Q, pol_llstat, HWE_llstat);
								for (ig=1;ig<nsample;ig++) {
									fprintf(outstream, "%s\t", quartet[ig].c_str());
									// printf("%s\t", quartet[ig].c_str());
								}
								fprintf(outstream, "%s\n", quartet[nsample].c_str());
								// printf("%s\n", quartet[nsample].c_str());
							} else if (s_mode == "F") {
								if (ref_info == 1) {
									if (best_error <= max_e) {
										fprintf(outstream, "%s\t%d\t%s\t%s\tNA\t%d\t%f\t%f\t%f\t%f\t%f\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), pop_coverage, Ni, best_p, best_q, best_error, best_H);
										// printf("%s\t%d\t%s\t%s\tNA\t%d\t%f\t%f\t%f\t%f\t%f\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), pop_coverage, Ni, best_p, best_q, best_error, best_H);
									} else {
										fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
										// printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
									}
								} else if (ref_info == 0) {
									if (best_error <= max_e) {
										fprintf(outstream, "%s\tNA\t%d\t%f\t%f\t%f\t%f\t%f\tNA\n", major_allele.c_str(), pop_coverage, Ni, best_p, best_q, best_error, best_H);
 										// printf("%s\tNA\t%d\t%f\t%f\t%f\t%f\t%f\tNA\n", major_allele.c_str(), pop_coverage, Ni, best_p, best_q, best_error, best_H);
									} else {
										fprintf(outstream, "NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
										// printf("NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
									}
								}
							} else if (s_mode == "B") {
								test_Mac = 2.0*nesample*best_p - (int)(2.0*nesample*best_p);
                                                                if ( test_Mac >= 0.5) {
                                                                        best_Mac = (int)(2.0*nesample*best_p) + 1;
                                                                } else {
                                                                        best_Mac = (int)(2.0*nesample*best_p);
                                                                }
                                                                test_mac = 2.0*nesample*best_q - (int)(2.0*nesample*best_q);
                                                                if ( test_mac >= 0.5) {
                                                                        best_mac = (int)(2.0*nesample*best_q) + 1;
                                                                } else {
                                                                        best_mac = (int)(2.0*nesample*best_q);
                                                                }
								if (ref_info == 1) {
                                                                	if (best_error <= max_e) {
                                                                        	fprintf(outstream, "%s\t%d\t%s\t%s\tNA\t%d\t%f\t%d\t%d\t%f\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), pop_coverage, Ni, best_Mac, best_mac, best_error);
                                                                                // printf("%s\t%d\t%s\t%s\tNA\t%d\t%f\t%d\t%d\t%f\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), pop_coverage, Ni, best_Mac, best_mac, best_error);
                                                                        } else {
                                                                        	fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
                                                                                // printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\n");
                                                                        }
                                                       		} else if (ref_info == 0) {
                                                                       if (best_error <= max_e) {
                                                                                fprintf(outstream, "%s\tNA\t%d\t%f\t%d\t%d\t%f\tNA\n", major_allele.c_str(), pop_coverage, Ni, best_Mac, best_mac, best_error);
                                                                                // printf("%s\tNA\t%d\t%f\t%d\t%d\t%f\tNA\n", major_allele.c_str(), pop_coverage, Ni, best_Mac, best_mac, best_error);
                                                                       } else {
                                                                                fprintf(outstream, "NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
                                                                                // printf("NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
                                                                       }
                                                                }
							} else if (s_mode == "p") {
								if (ref_info == 1) {
									if (best_error <= max_e) {
										fprintf(outstream, "%s\t%d\t%s\t%s\tNA\t%d\t%f\t%f\t%f\t%f\t%f\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), pop_coverage, Nc, best_p, best_q, best_error, best_H);
										// printf("%s\t%d\t%s\t%s\tNA\t%d\t%f\t%f\t%f\t%f\t%f\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), pop_coverage, Nc, best_p, best_q, best_error, best_H);
									} else {
										fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
										// printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
									}
								} else if (ref_info == 0) {
									if (best_error <= max_e) {
										fprintf(outstream, "%s\tNA\t%d\t%f\t%f\t%f\t%f\t%f\tNA\n", major_allele.c_str(), pop_coverage, Nc, best_p, best_q, best_error, best_H);
										// printf("%s\tNA\t%d\t%f\t%f\t%f\t%f\t%f\tNA\n", major_allele.c_str(), pop_coverage, Nc, best_p, best_q, best_error, best_H);
									} else {
										fprintf(outstream, "NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
										// printf("NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
									}
								}
							}
						}
					} else {	// ML estimates not found
						fprintf(stderr, "ML estimates not found at site %d on %s\n", site, scaffold.c_str());
						// Print out the results for f, c (with the as option equal to one), B, F, and p modes
						if (s_mode == "f") {
							fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
							// printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
						} else if (s_mode == "c" && as == 1) {
							fprintf(outstream, "%s\t%d\t%s\t%d\t%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t", scaffold.c_str(), site, ref_nuc.c_str(), n1, n2);
							// printf("%s\t%d\t%s\t%d\t%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t", scaffold.c_str(), site, ref_nuc.c_str(), n1, n2);
							for (ig=1;ig<nsample;ig++) {
								fprintf(outstream, "%s\t", quartet[ig].c_str());
								// printf("%s\t", quartet[ig].c_str());
							}
							fprintf(outstream, "%s\n", quartet[nsample].c_str());
							// printf("%s\n", quartet[nsample].c_str());
						} else if (s_mode == "F") {
							if (ref_info == 1) {
								fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
								// printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
							} else if (ref_info == 0) {
								fprintf(outstream, "NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
								// printf("NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
							}
						} else if (s_mode == "B") {
							if (ref_info == 1) {
								fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
								// printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\n");
							} else if (ref_info == 0) {
								fprintf(outstream, "NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
								// printf("NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
							}
						} else if (s_mode == "p") { 
							if (ref_info == 1) {
                                                                fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
                                                                // printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
                                                        } else if (ref_info == 0) {
                                                                fprintf(outstream, "NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
                                                                // printf("NA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\n", pop_coverage);
                                                        }
						}
					}	  
				}  	
			} else if (pop_read[n1] == pop_coverage) {
				// Print out the monomorphism results for the f, c (with the as option equal to one), B, F, and p modes
				if (s_mode == "f") {
					if (n1 == 1) {
						major_allele = "A";
					} else if (n1 == 2) {
        					major_allele = "C";
					} else if (n1 == 3) {
        					major_allele = "G";
					} else if (n1 == 4) {
        					major_allele = "T";
					}
					best_p = 1.0;
					best_q = 1.0 - best_p;
					test_Mac = 2.0*nesample*best_p - (int)(2.0*nesample*best_p);
					if ( test_Mac >= 0.5) {
                                		best_Mac = (int)(2.0*nesample*best_p) + 1;
                                	} else {
                                        	best_Mac = (int)(2.0*nesample*best_p);
                                	}
					test_mac = 2.0*nesample*best_q - (int)(2.0*nesample*best_q);
                                	if ( test_mac >= 0.5) {
                                		best_mac = (int)(2.0*nesample*best_q) + 1;
                                	} else {
                                        	best_mac = (int)(2.0*nesample*best_q);
                                	}
					best_error = 0.0;
					null_best_error = 0.0;
					best_D_A = 0.0;
					best_P = 1.0;
					best_H = 0.0;
					best_Q = 0.0;
					best_h = 0.0;
					fprintf(outstream, "%s\t%d\t%s\t%s\tNA\t%d\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\tNA\tNA\t%f\t%f\t%f\t%f\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), pop_coverage, Nc, Nr, best_Mac, best_mac, best_p, best_q, best_error, null_best_error, best_P, best_H, best_Q, best_h);
					// printf("%s\t%d\t%s\t%s\tNA\t%d\t%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\tNA\tNA\t%f\t%f\t%f\t%f\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), pop_coverage, Nc, Nr, best_Mac, best_mac, best_p, best_q, best_error, null_best_error, best_P, best_H, best_Q, best_h);
				} else if (s_mode == "c" && as == 1) {
					best_P = 1.0;
					best_H = 0.0;
					best_Q = 0.0;
					best_error = 0.0;
					fprintf(outstream, "%s\t%d\t%s\t%d\tNA\t%d\t%f\t%f\t%f\t%f\tNA\tNA\t", scaffold.c_str(), site, ref_nuc.c_str(), n1, pop_coverage, best_error, best_P, best_H, best_Q);
					// printf("%s\t%d\t%s\t%d\tNA\t%d\t%f\t%f\t%f\t%f\tNA\tNA\t", scaffold.c_str(), site, ref_nuc.c_str(), n1, pop_coverage, best_error, best_P, best_H, best_Q);
					for (ig=1;ig<nsample;ig++) {
						fprintf(outstream, "%s\t", quartet[ig].c_str());
						// printf("%s\t", quartet[ig].c_str());
					}
					fprintf(outstream, "%s\n", quartet[nsample].c_str());
					// printf("%s\n", quartet[nsample].c_str());	
				} else if (s_mode == "F") {
					best_p = 1.0;
					best_q = 0.0;
					best_error = 0.0;
                                        best_H = 0.0;
					if (ref_info == 1) {
						fprintf(outstream, "%s\t%d\t%s\t%s\tNA\t%d\t%f\t%f\t%f\t%f\t%f\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), pop_coverage, Ni, best_p, best_q, best_error, best_H);
						// printf("%s\t%d\t%s\t%s\tNA\t%d\t%f\t%f\t%f\t%f\t%f\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), pop_coverage, Ni, best_p, best_q, best_error, best_H);
					} else if (ref_info == 0) {
						fprintf(outstream, "%s\tNA\t%d\t%f\t%f\t%f\t%f\t%f\tNA\n", major_allele.c_str(), pop_coverage, Ni, best_p, best_q, best_error, best_H);
						// printf("%s\tNA\t%d\t%f\t%f\t%f\t%f\t%f\tNA\n", major_allele.c_str(), pop_coverage, Ni, best_p, best_q, best_error, best_H);
					}
				} else if (s_mode == "B") {
					if (n1 == 1) {
                                                major_allele = "A";
                                        } else if (n1 == 2) {
                                                major_allele = "C";
                                        } else if (n1 == 3) {
                                                major_allele = "G";
                                        } else if (n1 == 4) {
                                                major_allele = "T";
                                        }
                                        best_p = 1.0;
                                        best_q = 1.0 - best_p;
                                        test_Mac = 2.0*nesample*best_p - (int)(2.0*nesample*best_p);
                                        if ( test_Mac >= 0.5) {
                                                best_Mac = (int)(2.0*nesample*best_p) + 1;
                                        } else {
                                                best_Mac = (int)(2.0*nesample*best_p);
                                        }
                                        test_mac = 2.0*nesample*best_q - (int)(2.0*nesample*best_q);
                                        if ( test_mac >= 0.5) {
                                                best_mac = (int)(2.0*nesample*best_q) + 1;
                                        } else {
                                                best_mac = (int)(2.0*nesample*best_q);
                                        }
                                        best_error = 0.0;
					if (ref_info == 1) {
						fprintf(outstream, "%s\t%d\t%s\t%s\tNA\t%d\t%f\t%d\t%d\t%f\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), pop_coverage, Ni, best_Mac, best_mac, best_error);
						// printf("%s\t%d\t%s\t%s\tNA\t%d\t%f\t%d\t%d\t%f\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), pop_coverage, Ni, best_Mac, best_mac, best_error);
					} else if (ref_info == 0) {
						fprintf(outstream, "%s\tNA\t%d\t%f\t%d\t%d\t%f\tNA\n", major_allele.c_str(), pop_coverage, Ni, best_Mac, best_mac, best_error);
						// printf("%s\tNA\t%d\t%f\t%d\t%d\t%f\tNA\n", major_allele.c_str(), pop_coverage, Ni, best_Mac, best_mac, best_error);
					}
				} else if (s_mode == "p") {
                                        best_p = 1.0;
                                        best_q = 0.0;
                                        best_error = 0.0;
                                        best_H = 0.0;
                                        if (ref_info == 1) {
                                                fprintf(outstream, "%s\t%d\t%s\t%s\tNA\t%d\t%f\t%f\t%f\t%f\t%f\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), pop_coverage, Nc, best_p, best_q, best_error, best_H);
                                                // printf("%s\t%d\t%s\t%s\tNA\t%d\t%f\t%f\t%f\t%f\t%f\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), major_allele.c_str(), pop_coverage, Nc, best_p, best_q, best_error, best_H);
                                        } else if (ref_info == 0) {
                                                fprintf(outstream, "%s\tNA\t%d\t%f\t%f\t%f\t%f\t%f\tNA\n", major_allele.c_str(), pop_coverage, Nc, best_p, best_q, best_error, best_H);
                                                // printf("%s\tNA\t%d\t%f\t%f\t%f\t%f\t%f\tNA\n", major_allele.c_str(), pop_coverage, Nc, best_p, best_q, best_error, best_H);
                                        }
                                }		
			} 
		} else {
			// Print out the reference-sequence information when there is no data for the f, c (with the as option equal to one), B, F, and p modes
			if (s_mode == "f") {
				fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\t%f\t%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage, Nc, Nr);
				// printf("%s\t%d\t%s\tNA\tNA\t%d\t%f\t%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage, Nc, Nr);
			} else if (s_mode == "c" && as == 1) {
                                fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\t", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
                                // printf("%s\t%d\t%s\tNA\tNA\t%d\tNA\tNA\tNA\tNA\tNA\tNA\t", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage);
                                for (ig=1;ig<nsample;ig++) {
                                        fprintf(outstream, "%s\t", quartet[ig].c_str());
                                        // printf("%s\t", quartet[ig].c_str());
                                }
                                fprintf(outstream, "%s\n", quartet[nsample].c_str());
                                // printf("%s\n", quartet[nsample].c_str());
			} else if (s_mode == "F") {
				if (ref_info == 1) {
					fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\t%f\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage, Ni);
					// printf("%s\t%d\t%s\tNA\tNA\t%d\t%f\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage, Ni);
				} else if (ref_info == 0) {
					fprintf(outstream, "NA\tNA\t%d\t%f\tNA\tNA\tNA\tNA\tNA\n", pop_coverage, Ni);
					// printf("NA\tNA\t%d\t%f\tNA\tNA\tNA\tNA\tNA\n", pop_coverage, Ni);
				}
			} else if (s_mode == "B") {
				if (ref_info == 1) {
					fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\t%f\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage, Ni);
					// printf("%s\t%d\t%s\tNA\tNA\t%d\t%f\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage, Ni);
				} else if (ref_info == 0) {
					fprintf(outstream, "NA\tNA\t%d\t%f\tNA\tNA\tNA\tNA\n", pop_coverage, Ni);
					// printf("NA\tNA\t%d\t%f\tNA\tNA\tNA\tNA\n", pop_coverage, Ni);
				}
			} else if (s_mode == "p") {
                                if (ref_info == 1) {
                                        fprintf(outstream, "%s\t%d\t%s\tNA\tNA\t%d\t%f\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage, Nc);
                                        // printf("%s\t%d\t%s\tNA\tNA\t%d\t%f\tNA\tNA\tNA\tNA\tNA\n", scaffold.c_str(), site, ref_nuc.c_str(), pop_coverage, Nc);
                                } else if (ref_info == 0) {
                                        fprintf(outstream, "NA\tNA\t%d\t%f\tNA\tNA\tNA\tNA\tNA\n", pop_coverage, Nc);
                                        // printf("NA\tNA\t%d\t%f\tNA\tNA\tNA\tNA\tNA\n", pop_coverage, Nc);
                                }
                        }			
		}		
	}
	// Close the input and output files
	inputFile.close();
	fclose(outstream);
	
	return 0;
}
