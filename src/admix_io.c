#include <stdio.h>
#include <stdbool.h> 
#include <stdlib.h>
#include <sys/queue.h>
#include <getopt.h>
#include <dotgeno.h>

#define MAGIC_BYTES_SIZE 4
#define STR_BUF_EXTRA 6

bool HASH_CHECK = true;
bool IS_VERBOSE = false;

typedef enum {
    EGN,
	PAM
} geno_file_type;

typedef struct {
	char* snp;
	char* ind;
	char* geno;
} admixio_file_trio;

typedef union {
	egn_file_reader egn;
	pam_file_reader pam;
} geno_reader_base;

typedef union {
	egn_file_writer egn;
	pam_file_writer pam;
} geno_writer_base;

typedef struct {
	geno_reader_base reader;
	geno_file_type geno_type;
} geno_reader;

typedef struct {
	geno_writer_base writer;
	geno_file_type geno_type;
} geno_writer;

typedef struct {
	snp_data snp;
	ind_data ind;
	geno_reader geno;
} admixio_data_trio;

typedef struct {
	struct idx_head** elems;
	size_t length;
} idx_list_arr;

idx_list_arr init_idx_list_arr(size_t length) {
	idx_list_arr out;
	out.length = length;
	out.elems = (struct idx_head**)malloc(length * sizeof(struct idx_head*));
	return out;
}

void free_idx_list_arr(idx_list_arr* ila) {
	free(ila->elems);
}

geno_file_type get_geno_file_type(char* filename) {
	char magic_bytes[MAGIC_BYTES_SIZE + 1];
	magic_bytes[MAGIC_BYTES_SIZE] = '\0';
	FILE* fp = fopen(filename, "r");
	if(fp == NULL) {
		fprintf(stderr, "ERROR: cannot open file %s\n", filename);
		exit(EXIT_FAILURE);
	}
	size_t n_bytes_read = fread(magic_bytes, 1, MAGIC_BYTES_SIZE, fp);
	fclose(fp);
	// can only be a PAM if file size is greater than 4 bytes (MAGIC_BYTES_SIZE)
	if(n_bytes_read == MAGIC_BYTES_SIZE) {
		if(strcmp(magic_bytes, "GENO") == 0) {
			return PAM;
		}
	}
	// now check if it is an EGN if not a PAM
	bool is_egn = true;
	for(size_t i = 0; i < n_bytes_read; i++) {
		is_egn = is_egn && (magic_bytes[i] == '0' || magic_bytes[i] == '1' || magic_bytes[i] == '2' || magic_bytes[i] == '9');
	}
	if(is_egn) {
		return EGN;
	} else {
		fprintf(stderr, "ERROR: file %s is neither a PACKEDANCESTRYMAP nor an EIGENSTRAT file\n", filename);
		exit(EXIT_FAILURE);
	}
}

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

geno_writer writer_init(geno_file_type geno_type, char* filename, snp_data* snp_info, ind_data* ind_info) {
	geno_writer wtr;
	wtr.geno_type = geno_type;
	switch(geno_type) {
		case PAM:
			wtr.writer.pam = pam_file_writer_init(filename, snp_info, ind_info);
			write_pam_header(&wtr.writer.pam, snp_info, ind_info);
			break;
		case EGN:
			wtr.writer.egn = egn_file_writer_init(filename, snp_info, ind_info);
			break;
	}
	return wtr;
}

void write_record(geno_writer* wtr, uint8_t* dosages) {
	switch(wtr->geno_type) {
		case PAM:
			write_pam_record(&wtr->writer.pam, dosages);
			break;
		case EGN:
			write_egn_record(&wtr->writer.egn, dosages);
			break;
	}
}

uint8_t* read_record(geno_reader* rdr) {
	switch(rdr->geno_type) {
		case PAM:
			return read_pam_record(&rdr->reader.pam);
		case EGN:
			return read_egn_record(&rdr->reader.egn);
	}
}

short goto_var(geno_reader* rdr, snp_data* snp_info, char* var_name) {
	switch(rdr->geno_type) {
		case PAM:
			return goto_var_pam(&rdr->reader.pam, snp_info, var_name);
		case EGN:
			return goto_var_egn(&rdr->reader.egn, snp_info, var_name);
	}
}

void close_geno_reader(geno_reader* gr) {
	switch(gr->geno_type) {
		case PAM:
			close_pam_file_reader(&gr->reader.pam);
			break;
		case EGN:
			close_egn_file_reader(&gr->reader.egn);
			break;
	}
}

void close_geno_writer(geno_writer* gr) {
	switch(gr->geno_type) {
		case PAM:
			close_pam_file_writer(&gr->writer.pam);
			break;
		case EGN:
			close_egn_file_writer(&gr->writer.egn);
			break;
	}
}

size_t get_max_index(struct idx_node** arr, size_t length) {
        size_t max = 0;
        size_t max_idx = 0;
        for(size_t i = 0; i < length; i++) {
                if(arr[i]->idx >= max) {
                        max = arr[i]->idx;
                        max_idx = i;
                }
        }
        return max_idx;
}

bool all_equal(struct idx_node** arr, size_t length) {
        size_t first = arr[0]->idx;
        for(size_t i = 1; i < length; i++) {
                if(first != arr[i]->idx) {
                        return false;
                }
        }
        return true; 
}

/* now filter merge based on the linked lists */
// return length of output
// TEST!!!!
size_t intersect_idx(idx_list_arr* ila, struct idx_head* head_out) {
	size_t cnt = 0;
	// get first elements
	struct idx_node** cur_elems = (struct idx_node**)malloc(ila->length * sizeof(struct idx_node*));
	for(size_t i = 0; i < ila->length; i++) {
		cur_elems[i] = TAILQ_FIRST(ila->elems[i]);
	}
	
	// fill in output linked list with elements of first list
	struct idx_node* tmp_node;
	TAILQ_FOREACH(tmp_node, ila->elems[0], nodes) {
		struct idx_node* idn = (struct idx_node*)malloc(sizeof(struct idx_node));
		idn->idx = tmp_node->idx;
		TAILQ_INSERT_TAIL(head_out, idn, nodes);
	}
	struct idx_node* cur_elem_out = TAILQ_FIRST(head_out);

	while(true) {
		bool end_while = false;  // change value to exit while loop	
		if(all_equal(cur_elems, ila->length)) {
			cnt++;
			cur_elem_out = TAILQ_NEXT(cur_elem_out, nodes);
			for(size_t i = 0; i < ila->length; i++) {
				cur_elems[i] = TAILQ_NEXT(cur_elems[i], nodes);
				if((cur_elems[i] == NULL)) {
					if(i == 0) { end_while = true; break; }
					struct idx_node* old_val;
					while(cur_elem_out) {
						old_val = cur_elem_out;
						TAILQ_REMOVE(head_out, cur_elem_out, nodes);
						cur_elem_out = TAILQ_NEXT(cur_elem_out, nodes);
						free(old_val);
					}
					end_while = true;
					break;
				}
			}
		} else {
			size_t max_i = get_max_index(cur_elems, ila->length);
			for(size_t i = 0; i < ila->length; i++) {
				if(i == max_i) { continue; }
				if(cur_elems[i]->idx == cur_elems[max_i]->idx) { continue; }
				if(i == 0) {
					TAILQ_REMOVE(head_out, cur_elem_out, nodes);
					struct idx_node* old_val = cur_elem_out;
					cur_elem_out = TAILQ_NEXT(cur_elem_out, nodes);
					free(old_val);
				}
				cur_elems[i] = TAILQ_NEXT(cur_elems[i], nodes);
				if(cur_elems[i] == NULL) {
					if(i == 0) {
						end_while = true;
						break;
					} else {
						free_idx_list(head_out);
						end_while = true;
						break;
					}
				}
			}
		}
		if(end_while) {
			free(cur_elems);
			break;
		}
	}
	return cnt;
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

void print_help() {
	printf("Future help message here\n");
}

int main(int argc, char* argv[]) {
	static struct option long_options[] = {
		{"help",           no_argument,       NULL, 'h'},
		{"prefix",         required_argument, NULL, 'p'},
		{"geno",           required_argument, NULL, 'P'},
		{"snp",            required_argument, NULL, 's'},
		{"ind",            required_argument, NULL, 'i'},
		{"keep",           required_argument, NULL, 'k'},
		{"extract",        required_argument, NULL, 'e'},
		{"out",            required_argument, NULL, 'o'},
		{"output_type",    required_argument, NULL, 't'},
		{"ignore-hash",    no_argument,       0,    300},
		{"verbose",        no_argument,       0,    500},
		{0,                0,                 0,      0}
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
	geno_file_type output_type = PAM;
	while(1) {
		int c = getopt_long(argc, argv, "hp:P:s:i:k:e:o:", long_options, NULL);
		if(c == -1) { break; }
		switch(c) {
			case 'h':
				print_help();
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
				if(ind_filt_file) {
					fprintf(stderr, "ERROR: extract file of '%s' has already been specified!\n", ind_filt_file);
					exit(EXIT_FAILURE);
				}
				ind_filt_file = optarg;
				break;
			case 'k':
				if(snp_filt_file) {
					fprintf(stderr, "ERROR: keep file of '%s' has already been specified!\n", snp_filt_file);
					exit(EXIT_FAILURE);
				}
				snp_filt_file = optarg;
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
			case 300:
				HASH_CHECK = false;
				break;
			case 500:
				IS_VERBOSE = true;
				break;
			case '?':
				break;
			default:
				abort();
		}
	}
	// ADD HASH CHECK	
	admixio_data_trio adt = admixio_data_init(aft);
	if(adt.geno.geno_type == PAM) { read_pam_header(&adt.geno.reader.pam); }
	struct idx_head ind_idx_head;
	struct ind_idx_head missing_ind_idx;
	if(ind_filt_file) {
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
	// init writer
	// make outputs
	ind_data ind_out;
	snp_data snp_out;
	filter_ind_data(&adt.ind, &ind_out, &ind_idx_head);
	filter_snp_data(&adt.snp, &snp_out, &snp_idx_head);
	write_ind_data(&ind_out, out_files.ind);
	write_snp_data(&snp_out, out_files.snp);
	geno_writer wtr = writer_init(output_type, out_files.geno, &snp_out, &ind_out);
	struct idx_node* idx_snp;
	TAILQ_FOREACH(idx_snp, &snp_idx_head, nodes) {
		goto_var(&adt.geno, &adt.snp, adt.snp.var_id[idx_snp->idx]);
		uint8_t* dosages_in = read_record(&adt.geno);
		uint8_t* dosages_out = (uint8_t*)malloc(sizeof(uint8_t) * ind_out.length);
		size_t idx = 0;
		struct idx_node* idx_ind;
		TAILQ_FOREACH(idx_ind, &ind_idx_head, nodes) {
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
}
