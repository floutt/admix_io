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

/* hash table setup! */
/*! @brief comparison function for ind_idx hash. For use with khash.h */
#define ind_idx_equal(a, b) ((strcmp((a).ind_id, (b).ind_id) == 0) && strcmp((a).ind_pop, (b).ind_pop) == 0)

/* hash table stuff */
/**
 *  @brief Hash function for ind_idx for use with khash.h, based on the djb2 hash function
 *
 *  @param[in] ind_idx_struct ind_idx structure to be hashed
 *
 *  @return hash value of type khint_t
 */
static inline khint_t hash_ind_idx(ind_idx ind_idx_struct) {
        // String linker for .ind hashing
        int IND_LINK_LEN = 20;
        char IND_LINK[21] = "gzvrEy55bcEN0gqRqvL6";

        khint_t hash = 5381;
        
        // first string
        int len = strlen(ind_idx_struct.ind_id);
        for(int i = 0; i < len; i++) {
                hash = ((hash << 5) + hash) + ind_idx_struct.ind_id[i];
        }
        
        // linker string
        for(int i = 0; i < IND_LINK_LEN; i++) {
                hash = ((hash << 5) + hash) + IND_LINK[i];
        }

        // second string
        len = strlen(ind_idx_struct.ind_pop);
        for(int i = 0; i < len; i++) {
                hash = ((hash << 5) + hash) + ind_idx_struct.ind_pop[i];
        }
        return hash;
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

/* now filter merge based on the linked lists */
// BLAH
// BLAH
