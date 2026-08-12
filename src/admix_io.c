#include <stdio.h>
#include <stdbool.h> 
#include <stdlib.h>
#include <sys/queue.h>
#include <dotgeno.h>
#include "khash.h"

#define MAGIC_BYTES_SIZE 4

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

admixio_data_trio admixio_data_init(admixio_file_trio file_info) {
	admixio_data_trio out_adt;
	out_adt.snp = read_snp_file(file_info.snp);
	out_adt.ind = read_ind_file(file_info.ind);

	geno_reader rdr;
	rdr.geno_type = get_geno_file_type(file_info.geno);
	switch(rdr.geno_type) {
		case PAM:
			rdr.reader.pam = pam_file_reader_init(file_info.geno, &out_adt.snp, &out_adt.ind);
			break;
		case EGN:
			rdr.reader.egn = egn_file_reader_init(file_info.geno, &out_adt.snp, &out_adt.ind);
			break;
	}
	out_adt.geno = rdr;
	return out_adt;
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
		cur_elems[i] = STAILQ_FIRST(ila->elems[i]);
	}

	// fill in output linked list with elements of first list
	struct idx_node* tmp_node;
	STAILQ_FOREACH(tmp_node, ila->elems[0], nodes) {
		struct idx_node* idn = (struct idx_node*)malloc(sizeof(struct idx_node));
		idn->idx = tmp_node->idx;
		STAILQ_INSERT_TAIL(head_out, idn, nodes);
	}
	struct idx_node* cur_elem_out = STAILQ_FIRST(head_out);

	while(true) {
		bool end_while = false;  // change value to exit while loop	
		if(all_equal(cur_elems, ila->length)) {
			for(size_t i = 0; i < ila->length; i++) {
				cnt++;
				cur_elems[i] = STAILQ_NEXT(cur_elems[i], nodes);
				if(cur_elems[i] == NULL) {
					if(i == 0) { break; }
					struct idx_node* cur = cur_elem_out;
					struct idx_node* nxt;
					while(cur) {
						cur = STAILQ_NEXT(cur, nodes);
						if(cur == NULL) { break; }
						nxt = STAILQ_NEXT(cur, nodes);
						STAILQ_REMOVE(head_out, cur, idx_node, nodes);
						free(cur);
						cur = nxt;
					}
				}
			}
		} else {
			size_t max_i = get_max_index(cur_elems, ila->length);
			for(size_t i = 0; i < ila->length; i++) {
				if(i == max_i) { continue; }
				if(cur_elems[i]->idx == cur_elems[max_i]->idx) { continue; }
				if(i == 0) { STAILQ_REMOVE(head_out, cur_elem_out, idx_node, nodes); }
				cur_elems[i] = STAILQ_NEXT(cur_elems[i], nodes);
				if(cur_elems[i] == NULL) {
					free_idx_list(head_out);
					end_while = true;
					break;
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
