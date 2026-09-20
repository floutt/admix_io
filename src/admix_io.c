#include <stdio.h>
#include <stdbool.h> 
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <sys/queue.h>
#include <getopt.h>
#include <dotgeno.h>
#include "tailq_sort.h"
#include "io_generics.h"
#include "idx_funcs.h"

#define STR_BUF_EXTRA 6

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define CMP_IDX(a, b) ((a)->idx > (b)->idx)

bool HASH_CHECK = true;
bool IS_VERBOSE = false;

typedef struct {
	char* snp;
	char* ind;
	char* geno;
} admixio_file_trio;

typedef struct {
	snp_data snp;
	ind_data ind;
	geno_reader geno;
} admixio_data_trio;

admixio_data_trio admixio_data_init(admixio_file_trio aft) {
	admixio_data_trio out_adt;
	out_adt.snp = read_snp_file(aft.snp);
	out_adt.ind = read_ind_file(aft.ind);

	geno_reader rdr;
	rdr.geno_type = get_geno_file_type(aft.geno);
	switch(rdr.geno_type) {
		case PAM:
			rdr.reader.pam = pam_file_reader_init(aft.geno, &out_adt.snp, &out_adt.ind);
			break;
		case EGN:
			rdr.reader.egn = egn_file_reader_init(aft.geno, &out_adt.snp, &out_adt.ind);
			break;
	}
	out_adt.geno = rdr;
	return out_adt;
}

char** str_split(char* str, char delim, size_t* n_elems) {
	size_t str_len = strlen(str);
	size_t nc_buf = str_len > 0;
	// get num chars
	for(size_t i = 0; i < str_len; i++) {
		if(str[i] == delim) {
			nc_buf += 1;
		}
	}
	size_t* col_lens = (size_t*)malloc(nc_buf * sizeof(size_t)); 
	size_t cur_len = 0;
	size_t cur_col = 0;
	for(size_t i = 0; i < str_len; i++) {
		if(!(str[i] == delim)) {
			cur_len += 1;
		} else {
			col_lens[cur_col] = cur_len;
			cur_len = 0;
			cur_col += 1;
		}
	}
	col_lens[cur_col] = cur_len;
	char** out_arr = (char**)malloc(nc_buf * sizeof(char*));
	for(size_t i = 0; i < nc_buf; i++) {
		out_arr[i] = (char*)malloc((col_lens[i] + 1) * sizeof(char));
		out_arr[i][col_lens[i]] = '\0';
	}
	size_t subber = 0;
	cur_col = 0;
	for(size_t i = 0; i < str_len; i++) {
		if((str[i] == delim)) {
			cur_col += 1;
			subber = i + 1;
		} else {
			out_arr[cur_col][i - subber] = str[i];
		}
	}
	free(col_lens);
	*n_elems = nc_buf;
	return out_arr;
}

size_t num_lines(char* filename) {
	FILE* fp = fopen(filename, "r");
	if(fp == NULL) {
		fprintf(stderr, "ERROR: could not open file %s!\n", filename);
		exit(EXIT_FAILURE);
	}
	size_t n = 0;
	char ch;
	while ((ch = fgetc(fp)) != EOF) {
		if (ch == '\n') { n++; }
	}
	fclose(fp);
	return n;
}

void print_help(char* program) {
	printf("Usage: %s [OPTIONS]\n", program);
    printf("General options:\n");
    printf("\t-h, --help                                 Display help message and exit\n");
	printf("\t--ignore-hash                              Ignore hash check for PACKEDANCESTRYMAP file\n");
	printf("\t--verbose                                  Print verbose output\n");
    printf("Input options:\n");
    printf("\t-p, --prefix <prefix>                      Prefix of input files\n");
	printf("\t-g, --geno <filename>                      Input genotype file\n");
	printf("\t-s, --snp  <filename>                      Input SNP file\n");
	printf("\t-i, --ind <filename>                       Input individual file\n");
	printf("Output options:\n");
	printf("\t-o, --out <prefix>                         Output file prefix\n");
	printf("\t-t, --output-type {egn|pam}                Output file type\n");
	printf("Filter options:\n");
	printf("\t-k, --keep <filename>                      Keep individuals included in file. The file is a tab-separated file where the first and second column have individual and population IDs.\n");
	printf("\t-e, --extract <filename>                   Keep variants included in the file.\n");
	printf("\t--keep-pop <filename>                      Keep individuals in the populations specified in file\n");
	printf("\t-S, --sex <sex>                            Include individuals of specified sex\n");
	printf("\t-c, --chr <chrs>                           Comma-separated list of chrs to include\n");
	printf("\t-r, --range <ranges>                       Comma-separated list of position ranges to include. In the form {chr}:{start}-{end}\n");
	printf("\t--remove <filename>                        Remove individuals included in file. The file is a tab-separated file where the first and second column have individual and population IDs.\n");
	printf("\t--exclude <filename>                       Remove variants included in the file.\n");
	printf("\t--remove-pop <filename>                    Remove individuals in the populations specified in file\n");
	printf("\t--remove-sex <sex>                         Remove individuals of specified sex\n");
	printf("\t--remove-chr <chrs>                        Comma-separated list of chrs to remove\n");
	printf("\t--remove-range <ranges>                    Comma-separated list of position ranges to remove. In the form {chr}:{start}-{end}\n");
	printf("\t--maf {0..0.5}                             Keep variants with a minor allele frequency greater than or equal to specified value\n");
	printf("\t--max-maf {0..0.5}                         Keep variants with a minor allele frequency less than or equal to specified value\n");
	printf("\t--mac <int> (must be greater than 0)       Keep variants with a minor allele count greater than or equal to specified value\n");
	printf("\t--max-mac <int> (must be greater than 0)   Keep variants with a minor allele count less than or equal to specified value\n");
	printf("\t--msnp {0..1}                              Remove variants with a missingness rate greater than specified value\n");
}

double get_maf(uint8_t* dosages, struct idx_head* head) {
	uint64_t non_na_cnt = 0;
	uint64_t sum_dosage = 0;
	struct idx_node* idx;
	TAILQ_FOREACH(idx, head, nodes) {
		non_na_cnt += dosages[idx->idx] != NAN_VAL;
		sum_dosage += (dosages[idx->idx] != NAN_VAL) * dosages[idx->idx];
	}
	double af = sum_dosage / (2.0 * non_na_cnt);
	return MIN(af, fabs(1-af));
}

uint64_t get_mac(uint8_t* dosages, struct idx_head* head) {
	uint64_t sum_dosage = 0;
	uint64_t sum_total = 0;
	struct idx_node* idx;
	TAILQ_FOREACH(idx, head, nodes) {
		sum_dosage += (dosages[idx->idx] != NAN_VAL) * dosages[idx->idx];
		sum_total += (dosages[idx->idx] != NAN_VAL) * 2;
	}
	return MIN(sum_dosage, sum_total - sum_dosage);
}

double get_msnp(uint8_t* dosages, struct idx_head* head) {
	uint64_t na_cnt = 0;
	uint64_t total = 0;
	struct idx_node* idx;
	TAILQ_FOREACH(idx, head, nodes) {
		na_cnt += dosages[idx->idx] == NAN_VAL;
		total++;
	}
	return (1.0 * na_cnt) / total;
}

int main(int argc, char* argv[]) {
	static struct option long_options[] = {
		{"help",           no_argument,       NULL,  'h'},
		{"prefix",         required_argument, NULL,  'p'},
		{"geno",           required_argument, NULL,  'P'},
		{"snp",            required_argument, NULL,  's'},
		{"ind",            required_argument, NULL,  'i'},
		{"keep",           required_argument, NULL,  'k'},
		{"extract",        required_argument, NULL,  'e'},
		{"out",            required_argument, NULL,  'o'},
		{"output-type",    required_argument, NULL,  't'},
		{"sex",            required_argument, NULL,  'S'},
		{"chr",            required_argument, NULL,  'c'},
		{"range",          required_argument, NULL,  'r'},
		{"ignore-hash",    no_argument,       0,     300},
		{"keep-pop",       required_argument, NULL,  400},
		{"maf",            required_argument, NULL,  500},
		{"max-maf",        required_argument, NULL,  600},
		{"mac",            required_argument, NULL,  700},
		{"max-mac",        required_argument, NULL,  800},
		{"msnp",           required_argument, NULL,  900},
		{"verbose",        no_argument,       0,    1000},
		{"remove",         no_argument,       0,    1100},
		{"remove-pop",     required_argument, NULL, 1200},
		{"exclude",        required_argument, NULL, 1300},
		{"remove-sex",     required_argument, NULL, 1400},
		{"remove-chr",     required_argument, NULL, 1500},
		{"remove-range",   required_argument, NULL, 1600},
		{0,                0,                 0,       0}
	};
	
	admixio_file_trio aft;
	aft.ind = NULL;
	aft.snp = NULL;
	aft.geno = NULL;
	admixio_file_trio out_files;
	out_files.ind = NULL;
	out_files.snp = NULL;
	out_files.geno = NULL;
	char* ind_filt_file = NULL;
	char* snp_filt_file = NULL;
	char* population_file = NULL;
	char* sex = NULL;
	char* chr_str = NULL;
	char* range_str = NULL;
	char* ind_filt_file_neg = NULL;
	char* snp_filt_file_neg = NULL;
	char* population_file_neg = NULL;
	char* sex_neg = NULL;
	char* chr_str_neg = NULL;
	char* range_str_neg = NULL;
	double maf_min = 0;
	double maf_max = 0.5;
	double msnp = 0; // missingness rate
	uint64_t mac_min = 0;
	uint64_t mac_max = UINT64_MAX;
	geno_file_type output_type = PAM;
	size_t N_IND_FILT = 0;
	size_t N_SNP_FILT = 0;
	while(1) {
		int c = getopt_long(argc, argv, "hp:P:s:i:k:e:o:", long_options, NULL);
		if(c == -1) { break; }
		switch(c) {
			case 'h':
				print_help(argv[0]);
				exit(EXIT_SUCCESS);
			case 'p':
				if(aft.geno) {
					fprintf(stderr, "ERROR: geno files already provided!");
					exit(EXIT_FAILURE);
				}
				if(aft.snp) {
					fprintf(stderr, "ERROR: SNP files already provided!");
					exit(EXIT_FAILURE);
				}
				if(aft.ind) {
					fprintf(stderr, "ERROR: Individual files already provided!");
					exit(EXIT_FAILURE);
				}
				aft.snp = (char*)malloc(sizeof(char) * (strlen(optarg) + STR_BUF_EXTRA));
				aft.ind = (char*)malloc(sizeof(char) * (strlen(optarg) + STR_BUF_EXTRA));
				aft.geno = (char*)malloc(sizeof(char) * (strlen(optarg) + STR_BUF_EXTRA));
				sprintf(aft.snp, "%s.snp", optarg);
				sprintf(aft.ind, "%s.ind", optarg);
				sprintf(aft.geno, "%s.geno", optarg);
				break;
			case 'P':
				if(aft.geno) {
					fprintf(stderr, "ERROR: Multiple sets of geno files provided!\n");
					exit(EXIT_FAILURE);
				}
				aft.geno = strdup(optarg);
				break;
			case 's':
				if(aft.snp) {
					fprintf(stderr, "ERROR: Multiple sets of SNP files provided!\n");
					exit(EXIT_FAILURE);
				}
				aft.snp = strdup(optarg);
				break;
			case 'i':
				if(aft.ind) {
					fprintf(stderr, "ERROR: Multiple sets of individual files provided!\n");
					exit(EXIT_FAILURE);
				}
				aft.ind = strdup(optarg);
				break;
			case 'o':
				if(out_files.geno) {
					fprintf(stderr, "ERROR: output prefix has already been specified!\n");
					exit(EXIT_FAILURE);
				} else {
					out_files.ind = (char*)malloc(sizeof(char) * (strlen(optarg) + STR_BUF_EXTRA));
					out_files.snp = (char*)malloc(sizeof(char) * (strlen(optarg) + STR_BUF_EXTRA));
					out_files.geno = (char*)malloc(sizeof(char) * (strlen(optarg) + STR_BUF_EXTRA));
					sprintf(out_files.snp, "%s.snp", optarg);
					sprintf(out_files.ind, "%s.ind", optarg);
					sprintf(out_files.geno, "%s.geno", optarg);
				}
				break;
			case 'e':
				if(snp_filt_file) {
					fprintf(stderr, "ERROR: extract file of '%s' has already been specified!\n", snp_filt_file);
					exit(EXIT_FAILURE);
				}
				snp_filt_file = optarg;
				break;
			case 'k':
				if(ind_filt_file) {
					fprintf(stderr, "ERROR: keep file of '%s' has already been specified!\n", ind_filt_file);
					exit(EXIT_FAILURE);
				}
				ind_filt_file = optarg;
				break;
			case 't':
				if(strcmp(optarg, "pam") == 0) {
					output_type = PAM;
				} else if(strcmp(optarg, "egn") == 0) {
					output_type = EGN;
				} else {
					fprintf(stderr, "ERROR: invalid output type %s. Must be 'pam' or 'egn'.\n", optarg);
					exit(EXIT_FAILURE);
				}
				break;
			case 'S':
				if(sex) {
					fprintf(stderr, "ERROR: sex argument already provided!\n");
					exit(EXIT_FAILURE);
				}
				sex = optarg;
				break;
			case 'c':
				if(chr_str) {
					fprintf(stderr, "ERROR: chr argument already provided!\n");
					exit(EXIT_FAILURE);
				}
				chr_str = optarg;
				break;
			case 'r':
				if(range_str) {
					fprintf(stderr, "ERROR: chr argument already provided!\n");
					exit(EXIT_FAILURE);
				}
				range_str = optarg;
				break;
			case 300:
				HASH_CHECK = false;
				break;
			case 400:
				if(population_file) {
					fprintf(stderr, "ERROR: population file of '%s' has already been specified!\n",  population_file);
					exit(EXIT_FAILURE);
				}
				population_file = optarg;
				break;
			case 500:
				maf_min = strtod(optarg, NULL);
				if(!((maf_min >= 0) && (maf_min <= 0.5))) {
					fprintf(stderr, "ERROR: maf must be between 0 and 0.5!\n");
					exit(EXIT_FAILURE);
				}
				break;
			case 600:
				maf_max = strtod(optarg, NULL);
				if(!((maf_max >= 0) && (maf_max <= 0.5))) {
					fprintf(stderr, "ERROR: max-maf must be between 0 and 0.5!\n");
					exit(EXIT_FAILURE);
				}
				break;
			case 700:
				mac_min = strtoul(optarg, NULL, 10);
				break;
			case 800:
				mac_max = strtoul(optarg, NULL, 10);
				break;
			case 900:
				msnp = strtod(optarg, NULL);
				if(!((msnp > 0) && (msnp <= 1))) {
					fprintf(stderr, "ERROR: msnp must be greater than 0 and less than or equal to 1!\n");
					exit(EXIT_FAILURE);
				}
				break;
			case 1000:
				IS_VERBOSE = true;
				break;
			case 1100:
				if(ind_filt_file_neg) {
					fprintf(stderr, "ERROR: remove file of '%s' has already been specified!\n", ind_filt_file_neg);
					exit(EXIT_FAILURE);
				}
				ind_filt_file_neg = optarg;
				break;
			case 1200:
				if(population_file_neg) {
					fprintf(stderr, "ERROR: remove-pop file of '%s' has already been specified!\n",  population_file_neg);
					exit(EXIT_FAILURE);
				}
				population_file_neg = optarg;
				break;
			case 1300:
				if(snp_filt_file_neg) {
					fprintf(stderr, "ERROR: exclude file of '%s' has already been specified!\n", snp_filt_file_neg);
					exit(EXIT_FAILURE);
				}
				snp_filt_file_neg = optarg;
				break;
			case 1400:
				if(sex_neg) {
					fprintf(stderr, "ERROR: remove-sex argument already provided!\n");
					exit(EXIT_FAILURE);
				}
				sex_neg = optarg;
				break;
			case 1500:
				if(chr_str_neg) {
					fprintf(stderr, "ERROR: remove-chr argument already provided!\n");
					exit(EXIT_FAILURE);
				}
				chr_str_neg = optarg;
				break;
			case 1600:
				if(range_str_neg) {
					fprintf(stderr, "ERROR: remove_range argument already provided!\n");
					exit(EXIT_FAILURE);
				}
				range_str_neg = optarg;
				break;
			case '?':
				break;
			default:
				abort();
		}
	}

	if((aft.geno == NULL) || (aft.snp == NULL) || (aft.ind == NULL)) {
		fprintf(stderr, "ERROR: complete set of input files not provided!\n");
		exit(EXIT_FAILURE);
	}
	if((out_files.geno == NULL) || (out_files.snp == NULL) || (out_files.ind == NULL)) {
		fprintf(stderr, "ERROR: Output file prefix not provided!\n");
		exit(EXIT_FAILURE);
	}
	admixio_data_trio adt = admixio_data_init(aft);

	// hash check
	if(adt.geno.geno_type == PAM) {
		hdr_data hdr = read_pam_header(&adt.geno.reader.pam);
		if(HASH_CHECK && ((hdr.ind_hash != adt.ind.hash) || (hdr.snp_hash != adt.snp.hash))) {
			fprintf(stderr, "ERROR: Invalid hash in PACKEDANCESTRYMAP file!\n");
			exit(EXIT_FAILURE);	
		}
	}

	struct idx_head ind_idx_head;
	struct ind_idx_head missing_ind_idx;
	if(ind_filt_file) {
		N_IND_FILT++;
		TAILQ_INIT(&ind_idx_head);
		TAILQ_INIT(&missing_ind_idx);
		size_t n_elems = num_lines(ind_filt_file);
		char** ind_ids = (char**)malloc(sizeof(char*) * n_elems);
		char** ind_pops = (char**)malloc(sizeof(char*) * n_elems);
		FILE* fp = fopen(ind_filt_file, "r");
		char* line = NULL;
		size_t len = 0;
		ssize_t read;
		size_t i = 0;
		while ((read = getline(&line, &len, fp)) != -1) {
			line[read - 1] = '\0';
			size_t num_cols;
			char** elems = str_split(line, '\t', &num_cols);
			if(num_cols != 2) {
				fprintf(stderr, "ERROR: invalid number of columns in extract file %s. Expected two columns per line\n", ind_filt_file);
				exit(EXIT_FAILURE);
			}
			ind_ids[i] = strdup(elems[0]);
			ind_pops[i] = strdup(elems[1]);
			free(elems[0]);
			free(elems[1]);
			free(elems);
			i++;
		}
		free(line);
		fclose(fp);
		get_multiple_ind_idx(&adt.ind, ind_ids, ind_pops, n_elems, &ind_idx_head, &missing_ind_idx);
		if(IS_VERBOSE) {
			if(!TAILQ_EMPTY(&missing_ind_idx)) { 
				struct ind_idx_node* tmp_node;
				printf("The following individuals were not found in %s:\n", aft.ind);
				TAILQ_FOREACH(tmp_node, &missing_ind_idx, nodes) {
					printf("\tind: %s, pop: %s\n", tmp_node->iidx->ind_id, tmp_node->iidx->ind_pop);
				}
			}
		}
		for(size_t idx = 0; idx < n_elems; idx++) {
			free(ind_ids[idx]);
			free(ind_pops[idx]);
		}
		free(ind_ids);
		free(ind_pops);
	}

	struct idx_head snp_idx_head;
	struct str_list_head missing_snp_idx;
	if(snp_filt_file) {
		N_SNP_FILT++;
		TAILQ_INIT(&snp_idx_head);
		TAILQ_INIT(&missing_snp_idx);
		size_t n_elems = num_lines(snp_filt_file);
		char** snp_ids = (char**)malloc(sizeof(char*) * n_elems);
		FILE* fp = fopen(snp_filt_file, "r");
		char* line = NULL;
		size_t len = 0;
		ssize_t read;
		size_t i = 0;
		while ((read = getline(&line, &len, fp)) != -1) {
			line[read - 1] = '\0';
			snp_ids[i] = strdup(line);
			i++;
		}
		free(line);
		fclose(fp);
		get_multiple_snp_idx(&adt.snp, snp_ids, n_elems, &snp_idx_head, &missing_snp_idx);
		TAILQ_MERGESORT(&snp_idx_head, idx_head, idx_node, nodes, CMP_IDX);
		if(IS_VERBOSE) {
			if(!TAILQ_EMPTY(&missing_snp_idx)) { 
				struct str_node* tmp_node;
				printf("The following SNPs were not found in %s:\n", aft.snp);
				TAILQ_FOREACH(tmp_node, &missing_snp_idx, nodes) {
					printf("\t%s\n", tmp_node->str);
				}
			}
		}
		for(size_t idx = 0; idx < n_elems; idx++) {
			free(snp_ids[idx]);
		}
		free(snp_ids);
	}

	struct idx_head ind_pop_idx_head;
	struct str_list_head missing_ind_pops_idx;
	if(population_file) {
		N_IND_FILT++;
		TAILQ_INIT(&ind_pop_idx_head);
		size_t n_elems = num_lines(population_file);
		char** pops = (char**)malloc(sizeof(char*) * n_elems);
		FILE* fp = fopen(population_file, "r");
		char* line = NULL;
		size_t len = 0;
		ssize_t read;
		size_t i = 0;
		while ((read = getline(&line, &len, fp)) != -1) {
			line[read - 1] = '\0';
			pops[i] = strdup(line);
			i++;
		}
		free(line);
		fclose(fp);
		get_multiple_pops(&adt.ind, pops, n_elems, &ind_pop_idx_head, &missing_ind_pops_idx);
		if(IS_VERBOSE) {
			if(!TAILQ_EMPTY(&missing_ind_pops_idx)) { 
				struct str_node* tmp_node;
				printf("The following populations were not found in %s:\n", aft.snp);
				TAILQ_FOREACH(tmp_node, &missing_ind_pops_idx, nodes) {
					printf("\t%s\n", tmp_node->str);
				}
			}
		}
		for(size_t idx = 0; idx < n_elems; idx++) {
			free(pops[idx]);
		}
		free(pops);
	}

	struct idx_head ind_sex_idx_head;
	if(sex) {
		N_IND_FILT++;
		TAILQ_INIT(&ind_sex_idx_head);
		get_multiple_sex(&adt.ind, sex, &ind_sex_idx_head); 
	}

	struct idx_head snp_chr_idx_head;
	if(chr_str) {
		N_SNP_FILT++;
		TAILQ_INIT(&snp_chr_idx_head);
		size_t n_elems;
		char** chrs = str_split(chr_str, ',', &n_elems);
		get_multiple_chrs(&adt.snp, chrs, n_elems, &snp_chr_idx_head); 
		for(size_t i = 0; i < n_elems; i++) { free(chrs[i]); }
		free(chrs);
	}

	struct idx_head snp_range_idx_head;
	if(range_str) {
		N_SNP_FILT++;
		TAILQ_INIT(&snp_range_idx_head);
		size_t n_elems;
		char** ranges = str_split(range_str, ',', &n_elems);
		char** chrs = (char**)malloc(sizeof(char*) * n_elems);
		uint64_t* starts = (uint64_t*)malloc(sizeof(uint64_t) * n_elems);
		uint64_t* ends = (uint64_t*)malloc(sizeof(uint64_t) * n_elems);
		for(size_t i = 0; i < n_elems; i++) {
			size_t n_elems_prime;
			char* rng = ranges[i];
			char** chr_pos = str_split(rng, ':', &n_elems_prime);
			if(n_elems_prime != 2) {
				fprintf(stderr, "ERROR: invalid range string: '%s'.\n", rng);
				exit(EXIT_FAILURE);
			}
			char* chr = strdup(chr_pos[0]);
			char* start_end = chr_pos[1];
			char** start_end_split = str_split(start_end, '-', &n_elems_prime);
			if(n_elems_prime != 2) {
				fprintf(stderr, "ERROR: invalid range string '%s'.\n", rng);
				exit(EXIT_FAILURE);
			}
			uint64_t start = (uint64_t)strtoul(start_end_split[0], NULL, 10);
			uint64_t end = (uint64_t)strtoul(start_end_split[1], NULL, 10);
			chrs[i] = chr;
			starts[i] = start;
			ends[i] = end;
			// freeing time!
			free(start_end_split[0]); free(start_end_split[1]); free(start_end_split);
			free(chr_pos[0]); free(chr_pos[1]); free(chr_pos);
			free(ranges[i]);
		}
		free(ranges);
		get_multiple_ranges(&adt.snp, chrs, starts, ends, n_elems, &snp_range_idx_head);
		for(size_t i = 0; i < n_elems; i++) {
			free(chrs[i]);
		}
		free(chrs);
		free(starts);
		free(ends);
	}
		
	idx_list_arr snp_ila = init_idx_list_arr(N_SNP_FILT);
	size_t idx = 0;
	if(snp_filt_file) { snp_ila.elems[idx] = &snp_idx_head; idx++; }
	if(chr_str) { snp_ila.elems[idx] = &snp_chr_idx_head; idx++; }
	if(range_str) { snp_ila.elems[idx] = &snp_range_idx_head; idx++; }

	idx_list_arr ind_ila = init_idx_list_arr(N_IND_FILT);
	idx = 0;
	if(ind_filt_file) { ind_ila.elems[idx] = &ind_idx_head; idx++; }
	if(population_file) { ind_ila.elems[idx] = &ind_pop_idx_head; idx++; }
	if(sex) { ind_ila.elems[idx] = &ind_sex_idx_head; idx++; }

	struct idx_head ind_idx_final_head;
	struct idx_head snp_idx_final_head;
	TAILQ_INIT(&ind_idx_final_head);
	TAILQ_INIT(&snp_idx_final_head);

	if(N_IND_FILT > 0) {
		intersect_idx(&ind_ila, &ind_idx_final_head);
	}  else {
		insert_range(adt.ind.length, &ind_idx_final_head);
	}
	if(N_SNP_FILT > 0) {
		intersect_idx(&snp_ila, &snp_idx_final_head);
	} else {
		insert_range(adt.snp.length, &snp_idx_final_head);
	}
	
	// now do the negative fitlering stuff
	N_IND_FILT = 0;
	N_SNP_FILT = 0;

	struct idx_head ind_idx_head_neg;
	struct ind_idx_head missing_ind_idx_neg;
	if(ind_filt_file) {
		N_IND_FILT++;
		TAILQ_INIT(&ind_idx_head_neg);
		TAILQ_INIT(&missing_ind_idx_neg);
		size_t n_elems = num_lines(ind_filt_file_neg);
		char** ind_ids = (char**)malloc(sizeof(char*) * n_elems);
		char** ind_pops = (char**)malloc(sizeof(char*) * n_elems);
		FILE* fp = fopen(ind_filt_file, "r");
		char* line = NULL;
		size_t len = 0;
		ssize_t read;
		size_t i = 0;
		while ((read = getline(&line, &len, fp)) != -1) {
			line[read - 1] = '\0';
			size_t num_cols;
			char** elems = str_split(line, '\t', &num_cols);
			if(num_cols != 2) {
				fprintf(stderr, "ERROR: invalid number of columns in extract file %s. Expected two columns per line\n", ind_filt_file_neg);
				exit(EXIT_FAILURE);
			}
			ind_ids[i] = strdup(elems[0]);
			ind_pops[i] = strdup(elems[1]);
			free(elems[0]);
			free(elems[1]);
			free(elems);
			i++;
		}
		free(line);
		fclose(fp);
		get_multiple_ind_idx(&adt.ind, ind_ids, ind_pops, n_elems, &ind_idx_head_neg, &missing_ind_idx_neg);
		TAILQ_MERGESORT(&ind_idx_head, idx_head, idx_node, nodes, CMP_IDX);
		if(IS_VERBOSE) {
			if(!TAILQ_EMPTY(&missing_ind_idx_neg)) { 
				struct ind_idx_node* tmp_node;
				printf("The following individuals were not found in %s:\n", aft.ind);
				TAILQ_FOREACH(tmp_node, &missing_ind_idx_neg, nodes) {
					printf("\tind: %s, pop: %s\n", tmp_node->iidx->ind_id, tmp_node->iidx->ind_pop);
				}
			}
		}
		for(size_t idx = 0; idx < n_elems; idx++) {
			free(ind_ids[idx]);
			free(ind_pops[idx]);
		}
		free(ind_ids);
		free(ind_pops);
	}

	struct idx_head snp_idx_head_neg;
	struct str_list_head missing_snp_idx_neg;
	if(snp_filt_file_neg) {
		N_SNP_FILT++;
		TAILQ_INIT(&snp_idx_head_neg);
		TAILQ_INIT(&missing_snp_idx_neg);
		size_t n_elems = num_lines(snp_filt_file_neg);
		char** snp_ids = (char**)malloc(sizeof(char*) * n_elems);
		FILE* fp = fopen(snp_filt_file_neg, "r");
		char* line = NULL;
		size_t len = 0;
		ssize_t read;
		size_t i = 0;
		while ((read = getline(&line, &len, fp)) != -1) {
			line[read - 1] = '\0';
			snp_ids[i] = strdup(line);
			i++;
		}
		free(line);
		fclose(fp);
		get_multiple_snp_idx(&adt.snp, snp_ids, n_elems, &snp_idx_head_neg, &missing_snp_idx_neg);
		if(IS_VERBOSE) {
			if(!TAILQ_EMPTY(&missing_snp_idx_neg)) { 
				struct str_node* tmp_node;
				printf("The following SNPs were not found in %s:\n", aft.snp);
				TAILQ_FOREACH(tmp_node, &missing_snp_idx_neg, nodes) {
					printf("\t%s\n", tmp_node->str);
				}
			}
		}
		for(size_t idx = 0; idx < n_elems; idx++) {
			free(snp_ids[idx]);
		}
		free(snp_ids);
	}

	struct idx_head ind_pop_idx_head_neg;
	struct str_list_head missing_ind_pops_idx_neg;
	if(population_file_neg) {
		N_IND_FILT++;
		TAILQ_INIT(&ind_pop_idx_head_neg);
		size_t n_elems = num_lines(population_file_neg);
		char** pops = (char**)malloc(sizeof(char*) * n_elems);
		FILE* fp = fopen(population_file_neg, "r");
		char* line = NULL;
		size_t len = 0;
		ssize_t read;
		size_t i = 0;
		while ((read = getline(&line, &len, fp)) != -1) {
			line[read - 1] = '\0';
			pops[i] = strdup(line);
			i++;
		}
		free(line);
		fclose(fp);
		get_multiple_pops(&adt.ind, pops, n_elems, &ind_pop_idx_head_neg, &missing_ind_pops_idx_neg);
		if(IS_VERBOSE) {
			if(!TAILQ_EMPTY(&missing_ind_pops_idx_neg)) { 
				struct str_node* tmp_node;
				printf("The following populations were not found in %s:\n", aft.snp);
				TAILQ_FOREACH(tmp_node, &missing_ind_pops_idx_neg, nodes) {
					printf("\t%s\n", tmp_node->str);
				}
			}
		}
		for(size_t idx = 0; idx < n_elems; idx++) {
			free(pops[idx]);
		}
		free(pops);
	}

	struct idx_head ind_sex_idx_head_neg;
	if(sex_neg) {
		N_IND_FILT++;
		TAILQ_INIT(&ind_sex_idx_head_neg);
		get_multiple_sex(&adt.ind, sex_neg, &ind_sex_idx_head_neg); 
	}

	struct idx_head snp_chr_idx_head_neg;
	if(chr_str_neg) {
		N_SNP_FILT++;
		TAILQ_INIT(&snp_chr_idx_head_neg);
		size_t n_elems;
		char** chrs = str_split(chr_str_neg, ',', &n_elems);
		get_multiple_chrs(&adt.snp, chrs, n_elems, &snp_chr_idx_head_neg); 
		for(size_t i = 0; i < n_elems; i++) { free(chrs[i]); }
		free(chrs);
	}

	struct idx_head snp_range_idx_head_neg;
	if(range_str_neg) {
		N_SNP_FILT++;
		TAILQ_INIT(&snp_range_idx_head_neg);
		size_t n_elems;
		char** ranges = str_split(range_str_neg, ',', &n_elems);
		char** chrs = (char**)malloc(sizeof(char*) * n_elems);
		uint64_t* starts = (uint64_t*)malloc(sizeof(uint64_t) * n_elems);
		uint64_t* ends = (uint64_t*)malloc(sizeof(uint64_t) * n_elems);
		for(size_t i = 0; i < n_elems; i++) {
			size_t n_elems_prime;
			char* rng = ranges[i];
			char** chr_pos = str_split(rng, ':', &n_elems_prime);
			if(n_elems_prime != 2) {
				fprintf(stderr, "ERROR: invalid range string: '%s'.\n", rng);
				exit(EXIT_FAILURE);
			}
			char* chr = strdup(chr_pos[0]);
			char* start_end = chr_pos[1];
			char** start_end_split = str_split(start_end, '-', &n_elems_prime);
			if(n_elems_prime != 2) {
				fprintf(stderr, "ERROR: invalid range string '%s'.\n", rng);
				exit(EXIT_FAILURE);
			}
			uint64_t start = (uint64_t)strtoul(start_end_split[0], NULL, 10);
			uint64_t end = (uint64_t)strtoul(start_end_split[1], NULL, 10);
			chrs[i] = chr;
			starts[i] = start;
			ends[i] = end;
			// freeing time!
			free(start_end_split[0]); free(start_end_split[1]); free(start_end_split);
			free(chr_pos[0]); free(chr_pos[1]); free(chr_pos);
			free(ranges[i]);
		}
		free(ranges);
		get_multiple_ranges(&adt.snp, chrs, starts, ends, n_elems, &snp_range_idx_head_neg);
		for(size_t i = 0; i < n_elems; i++) {
			free(chrs[i]);
		}
		free(chrs);
		free(starts);
		free(ends);
	}

	idx_list_arr snp_ila_neg = init_idx_list_arr(N_SNP_FILT);
	idx = 0;
	if(snp_filt_file_neg) { snp_ila_neg.elems[idx] = &snp_idx_head_neg; idx++; }
	if(chr_str_neg) { snp_ila_neg.elems[idx] = &snp_chr_idx_head_neg; idx++; }
	if(range_str_neg) { snp_ila_neg.elems[idx] = &snp_range_idx_head_neg; idx++; }

	idx_list_arr ind_ila_neg = init_idx_list_arr(N_IND_FILT);
	idx = 0;
	if(ind_filt_file_neg) { ind_ila_neg.elems[idx] = &ind_idx_head_neg; idx++; }
	if(population_file_neg) { ind_ila_neg.elems[idx] = &ind_pop_idx_head_neg; idx++; }
	if(sex_neg) { ind_ila_neg.elems[idx] = &ind_sex_idx_head_neg; idx++; }

	if(N_IND_FILT > 0) { remove_idx(&ind_ila_neg, &ind_idx_final_head); }
	if(N_SNP_FILT > 0) { remove_idx(&snp_ila_neg, &snp_idx_final_head); }

	// check if empty
	if(TAILQ_EMPTY(&ind_idx_final_head)) {
		printf("WARNING: No individuals meet the filtering criteria specified. No output was produced.\n");
		return 0;
	}
	if(TAILQ_EMPTY(&snp_idx_final_head)) {
		printf("WARNING: No variants meet the filtering criteria specified. No output was produced.\n");
		return 0;
	}

	// snp filter walkthrough
	struct idx_node* idx_snp = TAILQ_FIRST(&snp_idx_final_head);
	while(1) {
		goto_var(&adt.geno, &adt.snp, adt.snp.var_id[idx_snp->idx]);
		uint8_t* dosages = read_record(&adt.geno);
		struct idx_node* idx_nxt = TAILQ_NEXT(idx_snp, nodes);
		struct idx_node* old_val = idx_snp;
		bool end_while = true;
		bool removed = false;
		if((maf_min > 0) || (maf_max < 0.5)) {
			end_while = false;
			double maf = get_maf(dosages, &ind_idx_final_head);
			if((maf < maf_min) || (maf > maf_max)) {
				TAILQ_REMOVE(&snp_idx_final_head, idx_snp, nodes);
				free(idx_snp);
				removed = true;
			}
		}

		if(!removed && ((mac_min != 0) || (mac_max != UINT64_MAX))) {
			end_while = false;
			uint64_t mac = get_mac(dosages, &ind_idx_final_head);
			if((mac < mac_min) || (mac > mac_max)) {
				TAILQ_REMOVE(&snp_idx_final_head, idx_snp, nodes);
				free(idx_snp);
				removed = true;
			}
		}
		
		if(!removed && (msnp > 0)) {
			end_while = false;
			double msnp_prime = get_msnp(dosages, &ind_idx_final_head);
			if(msnp_prime > msnp) {
				TAILQ_REMOVE(&snp_idx_final_head, idx_snp, nodes);
				free(idx_snp);
				removed = true;
			}
		}
		free(dosages);
		if((idx_nxt == NULL) || end_while) {
			break;
		}
		idx_snp = idx_nxt;
	}	

	if(TAILQ_EMPTY(&snp_idx_final_head)) {
		printf("WARNING: No variants meet the filtering criteria specified. No output was produced.\n");
		return 0;
	}

	// init writer
	// make outputs
	ind_data ind_out;
	snp_data snp_out;
	filter_ind_data(&adt.ind, &ind_out, &ind_idx_final_head);
	filter_snp_data(&adt.snp, &snp_out, &snp_idx_final_head);
	write_ind_data(&ind_out, out_files.ind);
	write_snp_data(&snp_out, out_files.snp);
	geno_writer wtr = writer_init(output_type, out_files.geno, &snp_out, &ind_out);
	// final write
	TAILQ_FOREACH(idx_snp, &snp_idx_final_head, nodes) {
		goto_var(&adt.geno, &adt.snp, adt.snp.var_id[idx_snp->idx]);
		uint8_t* dosages_in = read_record(&adt.geno);
		uint8_t* dosages_out = (uint8_t*)malloc(sizeof(uint8_t) * ind_out.length);
		size_t idx = 0;
		struct idx_node* idx_ind;
		TAILQ_FOREACH(idx_ind, &ind_idx_final_head, nodes) {
			dosages_out[idx] = dosages_in[idx_ind->idx];
			idx++;
		}
		write_record(&wtr, dosages_out);
		free(dosages_in);
		free(dosages_out);
	}
	close_geno_writer(&wtr);
	// free
	free(aft.ind);
	free(aft.snp);
	free(aft.geno);
	free(out_files.ind);
	free(out_files.snp);
	free(out_files.geno);
	free_snp_data(&adt.snp);
	free_snp_data(&snp_out);
	free_ind_data(&adt.ind);
	free_ind_data(&ind_out);
	close_geno_reader(&adt.geno);
	if(ind_filt_file) { free_idx_list(&ind_idx_head); }
	if(ind_filt_file) { free_ind_idx_list(&missing_ind_idx); }
	if(snp_filt_file) { free_idx_list(&snp_idx_head); }
	if(snp_filt_file) { free_str_list(&missing_snp_idx); }
	if(population_file) { free_idx_list(&ind_pop_idx_head); }
	if(population_file) { free_str_list(&missing_ind_pops_idx); }
	if(sex) { free_idx_list(&ind_sex_idx_head); }
	if(chr_str) { free_idx_list(&snp_chr_idx_head); }
	if(range_str) { free_idx_list(&snp_range_idx_head); }
	free_idx_list_arr(&snp_ila);
	free_idx_list_arr(&ind_ila);
	free_idx_list(&snp_idx_final_head);
	free_idx_list(&ind_idx_final_head);
	if(ind_filt_file_neg) { free_idx_list(&ind_idx_head_neg); }
	if(ind_filt_file_neg) { free_ind_idx_list(&missing_ind_idx_neg); }
	if(snp_filt_file_neg) { free_idx_list(&snp_idx_head_neg); }
	if(snp_filt_file_neg) { free_str_list(&missing_snp_idx_neg); }
	if(population_file_neg) { free_idx_list(&ind_pop_idx_head_neg); }
	if(population_file_neg) { free_str_list(&missing_ind_pops_idx_neg); }
	if(sex_neg) { free_idx_list(&ind_sex_idx_head_neg); }
	if(chr_str_neg) { free_idx_list(&snp_chr_idx_head_neg); }
	if(range_str_neg) { free_idx_list(&snp_range_idx_head_neg); }
	free_idx_list_arr(&snp_ila_neg);
	free_idx_list_arr(&ind_ila_neg);
	return 0;
}
