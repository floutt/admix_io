# admix_io
`admix_io` is a program for blazing-fast, [PLINK-style](https://www.cog-genomics.org/plink/2.0/filter) input filtering of PACKEDANCESTRYMAP and EIGENSTRAT files using [`dotgenolib`](https://github.com/floutt/dotgenolib/).

## Installation
The only dependency for `admix_io` is `dotgenolib`, instructions for installing that library can be found [here](https://github.com/floutt/dotgenolib/). `admix_io` can be installed by running:
```sh
make
```
and, to install it globally on your system, you can simply run (as superuser):
```sh
sudo make install
```
Local installation can be done by specifying the path to the directory containing the `bin` directory in which to put the binary. For instance, one who wishes to install the binary to `$HOME/bin` can simply run:
```sh
make install PREFIX=$HOME
```

## Usage
`admix_io` contains SNP and individual filtering options inspired by `plink`. For instance someone who wants to include variants in chromosome 1 which have a minor allele frequency greater than 10% in Mixtec and Nahua women can run the following command:
```sh
admix_io --prefix /path/to/large_genetics_dataset --chr 1 --maf 0.1 --sex F --keep-pop pops_file
```
Where `pops_file` contains the following:
```
Mixtec.DG
Nahua.DG
```
A more general overview of the available input filtering options can be found by running `admix_io -h`:
```
Usage: admix_io [OPTIONS]
General options:
	-h, --help                                 Display help message and exit
	--ignore-hash                              Ignore hash check for PACKEDANCESTRYMAP file
	--verbose                                  Print verbose output
Input options:
	-p, --prefix <prefix>                      Prefix of input files
	-g, --geno <filename>                      Input genotype file
	-s, --snp  <filename>                      Input SNP file
	-i, --ind <filename>                       Input individual file
Output options:
	-o, --out <prefix>                         Output file prefix
	-t, --output-type {egn|pam}                Output file type
Filter options:
	-k, --keep <filename>                      Keep individuals included in file. The file is a tab-separated file where the first and second column have individual and population IDs.
	-e, --extract <filename>                   Keep variants included in the file.
	--keep-pop <filename>                      Keep individuals in the populations specified in file
	-S, --sex <sex>                            Include individuals of specified sex
	-c, --chr <chrs>                           Comma-separated list of chrs to include
	-r, --range <ranges>                       Comma-separated list of position ranges to include. In the form {chr}:{start}-{end}
	--remove <filename>                        Remove individuals included in file. The file is a tab-separated file where the first and second column have individual and population IDs.
	--exclude <filename>                       Remove variants included in the file.
	--remove-pop <filename>                    Remove individuals in the populations specified in file
	--remove-sex <sex>                         Remove individuals of specified sex
	--remove-chr <chrs>                        Comma-separated list of chrs to remove
	--remove-range <ranges>                    Comma-separated list of position ranges to remove. In the form {chr}:{start}-{end}
	--maf {0..0.5}                             Keep variants with a minor allele frequency greater than or equal to specified value
	--max-maf {0..0.5}                         Keep variants with a minor allele frequency less than or equal to specified value
	--mac <int> (must be greater than 0)       Keep variants with a minor allele count greater than or equal to specified value
	--max-mac <int> (must be greater than 0)   Keep variants with a minor allele count less than or equal to specified value
	--msnp {0..1}                              Remove variants with a missingness rate greater than specified value
```

