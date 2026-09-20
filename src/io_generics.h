#include <dotgeno.h>

#define MAGIC_BYTES_SIZE 4

typedef enum {
    EGN,
	PAM
} geno_file_type;

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

