# Hybrid simulations from RADseq data

Scripts to identify and generate hybrids genotypes based on empirical parental alleles sequenced 
and processed using RADseq data. The scripts work with the default outputs of the Stacks software.

TODO: Update with additional usage info.

## Identify a set of "diagnostic" a parental SNPs

```sh
usage: make_diagnostic_snp_whitlist.py [-h] -s SUMSTATS -f FSTATS [-o OUTDIR] [-n N_SITES]
       [-e] [-a HWE_ALPHA] [-p] [-d MIN_FST] [-m MIN_MAF] [-r]

Identify a set of diagnostic SNPs between two parental populations based on the FSTATS and
SUMSTATS outputs of STACKS.

options:
  -h, --help            show this help message and exit
  -s SUMSTATS, --sumstats SUMSTATS
                        (str) Path to the SUMSTATS file from POPULATIONS.
  -f FSTATS, --fstats FSTATS
                        (str) Path to the FST_Y-Z file from POPULATIONS.
  -o OUTDIR, --outdir OUTDIR
                        Output directory. [str]
  -n N_SITES, --n-sites N_SITES
                        (int) Number of individuals simulated per population. [default=1000]
  -e, --hwe             (bool) Discard sites out of HWE, applied per population [default=False]
  -a HWE_ALPHA, --hwe-alpha HWE_ALPHA
                        (float) p-value cutoff to clasify a site as out of HWE [default=0.05]
  -p, --private         (bool) Keep only sites containing private alleles in the parental
                        populations [default=False]
  -d MIN_FST, --min-fst MIN_FST
                        (float) Min Fst value used to classify a site as divergent between
                        parental populations [default=0.25]
  -m MIN_MAF, --min-maf MIN_MAF
                        (float) Minimum allele frequency to retain a parental allele [default=0.05]
  -r, --write-random-snp
                        (bool) Export only one random SNP per locus [default=False]
```

## Simulate hybrid genotypes from pool of parental alleles

```sh
usage: sim_hybrids_from_parents_vcf.py [-h] -p POPMAP -v VCF [-o OUTDIR] [-g GENERATIONS]
                                       [-n N_INDIVIDUALS] [-m MIN_MAF]

Simulate sequential hybrid populations (e.g., F1s, F2 backcrosses, etc.) from a pool of
parental alleles from two populations/species.

options:
  -h, --help            show this help message and exit
  -p POPMAP, --popmap POPMAP
                        (str) Path to population map file.
  -v VCF, --vcf VCF     Path to parental alleles in VCF format.
  -o OUTDIR, --outdir OUTDIR
                        Output directory. [str]
  -g GENERATIONS, --generations GENERATIONS
                        (int) Number of hybrid generations to simulate. [default=10]
  -n N_INDIVIDUALS, --n-individuals N_INDIVIDUALS
                        (int) Number of individuals simulated per population. [default=10]
  -m MIN_MAF, --min-maf MIN_MAF
                        (float) Minimum allele frequency to retain a parental allele
                        [default=0.05]
```

## Author
Angel G. Rivera-Colon  
arcolon14@gmail.com
