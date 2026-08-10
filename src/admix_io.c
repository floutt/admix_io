#include <stdio.h>
#include <stdbool.h> 
#include <stdlib.h>
#include <dotgeno.h>

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

/* individual filter functions */
// BLAH
// BLAH

/* snp filter functions */
// BLAH
// BLAH

/* merge multiple filters for ind and snp respectively */
// BLAH
// BLAH

/* now filter based on the linked lists */
// BLAH
// BLAH
