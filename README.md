# rad_hybrid_sims
Simulate hybrid genotypes from a pool of empirical parental alleles.

TODO: Update with additional usage info.

```sh
usage: sim_hybrids_from_parents_vcf.py [-h] -p POPMAP -v VCF [-o OUTDIR] [-g GENERATIONS] [-n N_INDIVIDUALS]

Simulate sequential hybrid populations (e.g., F1s, F2 backcrosses, etc.) from a pool of parental alleles from two populations/species.

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
                        Number of individuals simulated per population. [int, default=10]
```

## Author
Angel G. Rivera-Colon
arcolon14@gmail.com
