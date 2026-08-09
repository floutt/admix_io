#include <stdio.h>
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
	egn_file_reader egn:
	pam_file_reader pam;
} geno_reader_ptr;

typedef union {
	egn_file_writer egn:
	pam_file_writer pam;
} geno_writer_ptr;

typedef struct {
	snp_data* snp:
	ind_data* ind:
	geno_reader_ptr rdr;
	geno_file_type geno_type;
}
