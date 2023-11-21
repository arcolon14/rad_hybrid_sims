# Hybrid simulations from RADseq data

Scripts to identify and generate hybrids genotypes based on empirical parental alleles sequenced 
and processed using RADseq data. The scripts work with the default outputs of the Stacks software.

TODO: Update with additional usage info.

## Simulate hybrid genotypes from pool of parental alleles

```sh
usage: sim_hybrids_from_parents_vcf.py [-h] -p POPMAP -v VCF [-o OUTDIR] [-g GENERATIONS]
                                       [-n N_INDIVIDUALS] [-m MIN_MAF]

Simulate sequential hybrid populations (e.g., F1s, F2 backcrosses, etc.) from a pool of 
parental alleles from two populations/species.

options:
  -h, --help            show this help message and exit
  -p POPMAP, --popmap POPMAP
                        Path to population map file. [str]
  -v VCF, --vcf VCF     Path to parental alleles in VCF format. [str]
  -o OUTDIR, --outdir OUTDIR
                        Output directory. [str]
  -g GENERATIONS, --generations GENERATIONS
                        Number of hybrid generations to simulate. [int, default=10]
  -n N_INDIVIDUALS, --n-individuals N_INDIVIDUALS
                        Number of individuals simulated per population. [int,
                        default=10]
  -m MIN_MAF, --min-maf MIN_MAF
                        Minimum allele frequency to retain a parental allele [float,
                        default=0.05]
```

## Author
Angel G. Rivera-Colon
arcolon14@gmail.com
