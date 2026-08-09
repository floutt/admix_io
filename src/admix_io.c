#include <stdio.h>
#include <stdbool.h> 
#include <dotgeno.h>

typedef enum {
    EGN,
	PAM
} geno_file_type;

typedef struct {
	char* snp:
	char* ind:
	char* geno:
	geno_file_type geno_type;
} admixio_file_trio;

typedef union {
	egn_file_reader* egn:
	pam_file_reader* pam;
} geno_reader_base;

typedef union {
	egn_file_writer* egn:
	pam_file_writer* pam;
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
	snp_data* snp:
	ind_data* ind:
	geno_reader geno;
} admixio_data_trio;

geno_file_type get_geno_file_type(char* filename) {
	char magic_bytes[4];
	FILE* fp = fopen(filename, "r");
	if(fp == NULL) {
		fprintf(stderr, "ERROR: cannot open file %s\n", filename);
		exit(EXIT_FAILURE);
	}
	size_t n_bytes_read = fread(magic_bytes, 1, 4, fp);
	// can only be a PAM if file size is greater than 4 bytes
	if(n_bytes_read == 4) {
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
	}
}

admixio_data_trio admixio_data_init(admixio_file_trio file_info) {
/// BLAH BLAH
}
