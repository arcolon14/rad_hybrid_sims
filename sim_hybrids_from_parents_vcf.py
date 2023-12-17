#!/usr/bin/env python3
import sys, os, argparse, gzip, random, datetime, math

# (c) 2023 Angel G. Rivera-Colon

PROG = sys.argv[0].split('/')[-1]
DESC = '''Simulate sequential hybrid populations (e.g., F1s, F2
backcrosses, etc.) from a pool of parental alleles from two
populations/species.'''

def parse_args():
    p = argparse.ArgumentParser(prog=PROG, description=DESC)
    p.add_argument('-p', '--popmap', required=True,
                   help='(str) Path to population map file.')
    p.add_argument('-v', '--vcf', required=True,
                   help='Path to parental alleles in VCF format.')
    p.add_argument('-o', '--outdir', required=False, default='.',
                   help='Output directory.  [str]')
    p.add_argument('-g', '--generations', required=False, default=10, type=int,
                   help='(int) Number of hybrid generations to simulate. [default=10]')
    p.add_argument('-n', '--n-individuals', required=False, default=10, type=int,
                   help='(int) Number of individuals simulated per population. [default=10]')
    p.add_argument('-m', '--min-maf', required=False, default=0.05, type=float,
                   help='(float) Minimum allele frequency to retain a parental allele [default=0.05]')
    # Check input arguments
    args = p.parse_args()
    args.outdir = args.outdir.rstrip('/')
    if not os.path.exists(args.outdir):
        sys.exit(f"Error: '{args.outdir}' not found.")
    if not os.path.exists(args.popmap):
        sys.exit(f"Error: '{args.popmap}' not found.")
    if not os.path.exists(args.vcf):
        sys.exit(f"Error: '{args.vcf}' not found.")
    if args.generations <= 0:
        sys.exit(f"Error: generations ({args.generations}) must be greater than 0.")
    if args.n_individuals <= 0:
        sys.exit(f"Error: number of individuals ({args.n_individuals}) must be greater than 0.")
    if not 0.0 <= args.min_maf <= 0.5:
        sys.exit(f"Error: Min minor allele frequency ({args.min_maf}) must be between 0.0 - 0.5.")
    return args

#
# Classes
#

class SnpInfo():
    def __init__(self, chrom, bp, snpid, ref, alt):
        self.chr = chrom
        self.bp  = bp
        self.id  = snpid
        self.ref = ref
        self.alt = alt
    def __str__(self):
        return f'{self.chr} {self.bp} {self.id} {self.ref} {self.alt}'

class Cross():
    def __init__(self, generation, parent1, parent2, par1_gen=0, par2_gen=0):
        assert type(generation) is int
        assert generation >= 0
        assert type(par1_gen) is int
        if generation > 0:
            assert par1_gen == generation-1
        assert type(par2_gen) is int
        if generation > 0:
            assert 0 <= par2_gen <= par1_gen
        parentals = {'P1','P2'}
        self.gen       = generation
        self.par1      = parent1
        self.par1_g    = par1_gen # Generation of parent 1
        self.par2      = parent2
        self.par2_g    = par2_gen # Generation of parent 2
        self.hybrid    = True # Hybrid or not
        self.pop       = f'F{self.gen}'
        self.backcross = False # Backcross or not
        self.bc_dir    = None # Direction of backcross
        # Not a hybrid if the cross is between parentals
        if (self.par1 in {'P1','P2'}) and (self.par1 == self.par2):
            self.hybrid = False
            self.pop  = self.par1
            assert self.gen == 0, 'Error: Parental crosses must be generation 0'
            assert self.par1_g == 0, 'Error: Parental crosses must be generation 0'
            assert self.par2_g == 0, 'Error: Parental crosses must be generation 0'
        # Only pure parental can be gen 0
        if self.gen == 0:
            assert self.hybrid is False, 'Error: Generation 0 must be pure parentals'
        # F1s must only be from P1xP2
        elif self.gen == 1:
            F1 = True
            if self.par1 not in parentals:
                F1 = False
            elif self.par2 not in parentals:
                F1 = False
            elif self.par1 == self.par2:
                F1 = False
            assert F1, f'Error: For F1s but be a P1xP2 cross ({self.par1} x {self.par2}).'
        # Process later-generation crosses.
        else:
            # Late gen crosses cannot be between parentals
            if self.par1 in parentals and self.par2 in parentals:
                    assert False, f'For generations >1 both parents cannot both be parents ({self.par1}, {self.par2})'
            # For crosses among groups
            if self.par1 != self.par2:
                # Check for backcrosses to par1
                if self.par1 in parentals:
                    self.bc_dir = self.par1
                    self.pop = f'F{self.gen}-{self.bc_dir}'
                    self.backcross = True
                # Check for backcrosses to par2
                elif self.par2 in parentals:
                    self.bc_dir = self.par2
                    self.pop = f'F{self.gen}-{self.bc_dir}'
                    self.backcross = True
                # Crosses among hybrid groups
                else:
                    self.pop = f'F{self.gen}-H{self.par1_g}.{self.par2_g}'
            # For other Hybrid x Hybrid crosses
            else:
                self.pop = f'F{self.gen}-H{self.par1_g}.{self.par2_g}'
    def __str__(self):
        return f'{self.pop} {self.par1}x{self.par2} {self.hybrid} {self.backcross} {self.bc_dir} {self.par1_g}x{self.par2_g}'

class SampleAncestry:
    def __init__(self, sample_id, sample_pop, parent1, parent2, generation):
        assert parent1 != parent2, 'Error: parents cannot be the same individual.'
        self.id   = sample_id
        self.pop  = sample_pop
        self.par1 = parent1
        self.par2 = parent2
        self.gen  = generation
    def __str__(self):
        return f'{self.id} {self.pop} {self.par1} {self.par2} {self.gen}'

def now():
    '''Print the current date and time.'''
    return f'{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'

def time():
    '''Print today's date'''
    return f'{datetime.datetime.now().strftime("%Y%m%d")}'

def initialize_log(args, log=sys.stdout):
    '''Initialize the log file'''
    print(f'''{PROG} started on {now()}\n
Input parameters:
    Popmap:        {args.popmap}
    VCF:           {args.vcf}
    Outdir:        {args.outdir}
    Generations:   {args.generations}
    N individuals: {args.n_individuals}
    Min MAF:       {args.min_maf:0.06g}''',
    file=log, flush=True)

def parse_popmap(popmap_f, log=sys.stdout):
    '''Parse the popmap file and split individuals across parental populations.'''
    print('\nParsing parental popmap...', file=log, flush=True)
    parents = dict()
    n_sams = dict()
    with open(popmap_f) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.strip('\n').split('\t')
            sample = fields[0]
            pop = fields[1]
            if pop not in {'P1', 'P2'}:
                sys.exit("Error: population IDs in popmap must be 'P1' (parent 1) and 'P2' (parent 2).")
            parents[sample] = pop
            n_sams.setdefault(pop, 0)
            n_sams[pop] += 1
    # Report
    print(f'    Read {sum(n_sams.values()):,} individuals from the popmap file across {len(n_sams.keys())} populations:', file=log)
    for pop in sorted(n_sams):
        print(f'    {pop}: {n_sams[pop]:,} individuals', file=log)
    return parents

def parse_vcf(vcf_f, parents, vcf_out, min_maf=0.05, log=sys.stdout):
    '''Parse the VCF to obtain the parental allele counts/freqs.'''
    print('\nParsing VCF...', file=log, flush=True)
    n_snps   = 0
    vcf_rec  = 0
    parental_allele_freqs = dict()
    snp_info = dict()
    vcf_samples = None
    # Parse VCF
    with gzip.open(vcf_f, 'rt') if vcf_f.endswith('.gz') else open(vcf_f) as fh:
        for line in fh:
            # The header lines need to be added to the new output VCF file
            if line.startswith('##'):
                if line.startswith('##fileDate'):
                    vcf_out.write(f'##fileDate={time()}\n')
                elif line.startswith('##source'):
                    vcf_out.write(f'##source="{PROG}"\n')
                elif line.startswith('##INFO'):
                    continue
                elif line.startswith('##FORMAT'):
                    continue
                else:
                    vcf_out.write(line)
                continue
            elif line.startswith('#CHROM'):
                vcf_samples = line.strip('\n').split('\t')[9:]
                if len(vcf_samples) != len(parents):
                    sys.exit(f"Error: Number of samples in popmap ({len(parents)}) and input VCF ({len(vcf_samples)}) do not match.")
                print(f'    Processing genotypes from {len(vcf_samples):,} parental samples.', file=log, flush=True)
                continue
            else:
                # Process all other rows which contain variant sites
                vcf_rec += 1
                # Extract the reference info for each variant site,
                # which is used to reconstruct the final VCF
                fields = line.strip('\n').split('\t')
                chrom = fields[0]
                basepair = fields[1]
                id = fields[2]
                ref = fields[3]
                alt = fields[4]
                info = SnpInfo(chrom, basepair, id, ref, alt)
                # Process the genotype columns
                snp_allele_cnt = dict()
                for i, geno in enumerate(fields[9:]):
                    sample = vcf_samples[i]
                    pop = parents.get(sample, None)
                    if pop is None:
                        sys.exit(f"Error: sample VCF '{sample}' not found in popmap.")
                    geno = geno.split(':')[0]
                    alleles = geno.split('/')
                    # Initialize the allele count dictionary
                    # Tally the alleles seen in each population:
                    # { P1 : [ ref_cnt, alt_cnt ], P2 : [ ref_cnt, alt_cnt ] }
                    snp_allele_cnt.setdefault(pop, [0, 0])
                    for allele in alleles:
                        if allele.isnumeric():
                            allele = int(allele)
                            snp_allele_cnt[pop][allele] += 1
                # Convert the allele counts to frequencies
                snp_allele_freq = dict()
                for pop in snp_allele_cnt:
                    pop_cnts = snp_allele_cnt[pop]
                    pop_total = sum(pop_cnts)
                    assert pop_total > 0, f"Error: allele count for population {pop} must be greater than 0"
                    pop_freq = [ cnt/pop_total for cnt in pop_cnts ]
                    # Filter if under a MAF cutoff
                    minor_allele = min(pop_freq)
                    if minor_allele == 0.0 or minor_allele >= min_maf:
                        snp_allele_freq[pop] = pop_freq
                # Add to the genome-wide dictionary
                # Only when the two parents are seen
                if len(snp_allele_freq) == 2:
                    parental_allele_freqs[info.id] = snp_allele_freq
                    snp_info[id] = info
                n_snps += 1
        # Report and Return
        print(f'''    Read {vcf_rec:,} records from input VCF.
    Retained {len(snp_info):,} variant sites after applying filters.''', file=log, flush=True)
        return parental_allele_freqs, snp_info

def generate_crosses(generations, log=sys.stdout, ngen_hyb_hyb=4):
    '''Generate and check the crosses of interest. Return a list of crosses.'''
    assert generations > 0
    assert ngen_hyb_hyb <= generations
    crosses = list()
    # Make a list of possible hybrid parents for hybrid x hybrid crosses
    hybrid_parents = list()
    for g in range(0, generations+1):
        # Parental crosses
        if g == 0:
            for par in ['P1', 'P2']:
                cross = Cross(g, par, par, g, g)
                crosses.append(cross)
        # F1s
        elif g == 1:
            cross = Cross(g, 'P1', 'P2', 0, 0)
            crosses.append(cross)
            hybrid_parents.append(f'F{g-1}')
        # F2s
        elif g == 2:
            par1 = 'F1' # Parent 1 must always be an F1
            # F1 x P1 cross:
            cross = Cross(g, par1, 'P1', g-1, 0)
            crosses.append(cross)
            # F1 x P2 cross:
            cross = Cross(g, par1, 'P2', g-1, 0)
            crosses.append(cross)
            # F1 x F1 cross
            cross = Cross(g, par1, 'F1', g-1, g-1)
            crosses.append(cross)
        # Other hybrids, Fg
        else:
            # For a hybrid of generation g, one parent
            # MUST be from the previous generation.
            par1 = f'F{g-1}'
            # The other parent can be from any prior generation.
            # We call this k, ranging from 0 to g-1.
            for j in range(g):
                # Backcrosses (Fg x Parent)
                if j == 0:
                    for par in ['P1', 'P2']:
                        # These are "perfect" backcrosses,
                        # always against the same parent
                        par1 = f'F{g-1}-{par}'
                        cross = Cross(g, par1, par, g-1, j)
                        crosses.append(cross)
                # Crosses against the same generation hybrid
                # Fg x Fg
                elif j == (g-1):
                    par = f'F{g-1}-H{g-2}.{g-2}'
                    cross = Cross(g, par, par, g-1, g-1)
                    crosses.append(cross)
                # For all other crosses:
                else:
                    # Only process secondary hybrid x hybrid
                    # crosses for a certain number of generations
                    # else it creates excessive combinations.
                    if g > ngen_hyb_hyb:
                        continue
                    # Make crosses against F1s
                    # Fg x F1
                    if j == 1:
                        par1 = f'F{g-1}-H{g-2}.{g-2}'
                        par2 = 'F1'
                        cross = Cross(g, par1, par2, g-1, j)
                        crosses.append(cross)
                    # Process the other secondary hybrid crosses
                    # Fg x Fj
                    else:
                        par1 = f'F{g-1}-H{g-2}.{g-2}'
                        par1 = f'F{j}-H{j-1}.{j-1}'
                        cross = Cross(g, par1, par2, g-1, j)
                        crosses.append(cross)
    return crosses

def determine_parents(crosses, population_map):
    '''For a given map of individual and crosses, determine the ancestry of each
    specific sample (the id of its two parents).'''
    ancestry_map = dict()
    # The crosses MUST be sorted by generation
    for cross in sorted(crosses, key=lambda c: c.gen):
        assert isinstance(cross, Cross)
        ancestry_map.setdefault(cross.pop, [])
        generation = cross.gen
        pop_samples = population_map[cross.pop]
        # Loop over samples
        for sample in pop_samples:
            ancestry = None
            if generation == 0:
                # For parentals, their parents are from the same generation (you sample 2)
                parents = [sample, None]
                while sample in parents:
                    parents = random.sample(pop_samples, 2)
                ancestry = SampleAncestry(sample, cross.pop, parents[0], parents[1], generation)
                ancestry_map[cross.pop].append(ancestry)
            else:
                if cross.par1 != cross.par2:
                    # print(cross.par1, cross.par2)
                    # For all other pops, you sample for each population separately
                    parent1 = random.sample(population_map[cross.par1],1)[0]
                    parent2 = random.sample(population_map[cross.par2],1)[0]
                    ancestry = SampleAncestry(sample, cross.pop, parent1, parent2, generation)
                    ancestry_map[cross.pop].append(ancestry)
                else:
                    # This shouldn't happen without the hyb x hyb crosses, but just a check
                    parents = random.sample(population_map[cross.par1],2)
                    ancestry = SampleAncestry(sample, cross.pop, parents[0], parents[1], generation)
                    ancestry_map[cross.pop].append(ancestry)
    return ancestry_map

def make_vcf_header(crosses, ancestry_map, vcf):
    '''Re-make the VCF header including the sample IDs'''
    # Add some additional lines
    vcf.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
    vcf.write('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n')
    vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    header = ['#CHROM', 'POS', 'ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
    for cross in sorted(crosses, key=lambda c: c.gen):
        for sam in ancestry_map[cross.pop]:
            header.append(sam.id)
    header_str = '\t'.join(header)
    vcf.write(f'{header_str}\n')

def map_simulated_crosses(generations, n_individuals, outdir='.', log=sys.stdout, ngen_hyb_hyb=3):
    '''Make a map of the new hybrid populations alongside the IDs for simulated
    individuals. Write to disk as a popmap.tsv file. Return a list with the
    crosses and map of the ancestry of each individual (id of each parent).'''
    print('\nCreating a population map for simulated individuals...', file=log, flush=True)
    fh = open(f'{outdir}/simulated_hybrids.popmap.tsv', 'w')
    population_map = dict()
    # Calculate the possible crosses
    crosses = generate_crosses(generations, log)
    tot_pops = len(crosses)
    tot_sams = tot_pops * n_individuals
    print(f'    Simulating crosses for {tot_pops:,} populations.\n    A total of {tot_sams:,} individuals.', file=log, flush=True)
    # Loop over the crosses and print
    sample_i = 0
    for cross in crosses:
        assert isinstance(cross, Cross)
        # First, define the population ID:
        pop_id = cross.pop
        population_map.setdefault(pop_id, [])
        # Loop over the sampled individuals per pop
        for n in range(n_individuals):
            # The number of the individual in the population
            pop_sam_i = f'{n+1}'.zfill(len(str(n_individuals))) # Pad with Zeroes
            # Number of individual in the dataset
            sample_i += 1
            tot_sam_i = f'{sample_i}'.zfill(len(str(tot_sams)))
            # The sample ID combining population ID and sample number
            sample_id = f'{pop_id}_{pop_sam_i}_{tot_sam_i}'
            # Write to popmap
            fh.write(f'{sample_id}\t{pop_id}\n')
            population_map[pop_id].append(sample_id)
    fh.close()
    # Now, map the specific ancestry of each individual
    ancestry_map = determine_parents(crosses, population_map)
    return crosses, ancestry_map

def simulate_parental_site_genotypes(samples, popid, site_allele_freqs):
    '''Simulate the parental genotypes at a given site. Each genotype is
    sampled from the empirical parental allele frequency.'''
    pop_genotypes = dict()
    for sam in samples:
        assert isinstance(sam, SampleAncestry)
        pop_allele_freq = site_allele_freqs.get(popid, None)
        assert pop_allele_freq is not None, 'Error: parental allele freqs not found'
        genotype = random.choices([0,1], weights=pop_allele_freq, k=2)
        pop_genotypes[sam.id] = genotype
    return pop_genotypes

def sample_genotypes_from_allele_freqs(samples, cross, site_allele_freqs):
    '''For a given population, sample genotypes for a number of individuals
    according to the allele frequency observed in the population.'''
    assert isinstance(cross, Cross)
    pop_genotypes = dict()
    for sam in samples:
        assert isinstance(sam, SampleAncestry)
        # Sample allele from first parent
        al_p1 = sample_allele_from_freq(cross.par1, site_allele_freqs)
        # Sample allele from other parent
        al_p2 = sample_allele_from_freq(cross.par2, site_allele_freqs)
        pop_genotypes[sam.id] = [al_p1, al_p2]
    return pop_genotypes

def sample_allele_from_freq(popid, site_allele_freqs):
    '''Sample an allele from a population, given their allele freqs'''
    pop_allele_freq = site_allele_freqs.get(popid, None)
    assert pop_allele_freq is not None, 'Error: parental allele freqs not found'
    allele = random.choices([0,1], weights=pop_allele_freq, k=1)[0]
    return allele

def sample_hybrid_genotypes(samples, cross, site_genotypes):
    '''Sample genotypes for a hybrid population based on an existing pool 
    of discrete (existing) parental genotypes for a given site.'''
    pop_genotypes = dict()
    for sam in samples:
        assert isinstance(sam, SampleAncestry)
        # Get genotypes for parent 1
        par1_genotypes = site_genotypes[cross.par1][sam.par1]
        # Get genotypes for parent 2
        par2_genotypes = site_genotypes[cross.par2][sam.par2]
        # Sample parental alleles to get the offspring genotype
        genotype = [ random.sample(par1_genotypes,1)[0],
                     random.sample(par2_genotypes,1)[0] ]
        pop_genotypes[sam.id] = genotype
    return pop_genotypes

def format_vcf_row(site_info, site_genotypes, crosses, vcf):
    '''Aggregate the coordinate and allele information for each site, as well as the individual
    genotypes and format into the text of a VCF row. Save this now row text into the VCF fh.'''
    assert isinstance(site_info, SnpInfo)
    # Prepare the row
    row = [site_info.chr,  # CHROM
           site_info.bp,   # POS
           site_info.id,   # ID
           site_info.ref,  # REF
           site_info.alt,  # ALT
           '.',            # QUAL
           'PASS',         # FILTER
           'info',         # INFO (will be updated)
           'GT']           # FORMAT (we only have genotypes)
    # Two variables for the INFO field:
    info_ns  = 0  # Number of seen samples for INFO (should be all)
    info_af  = 0  # Allele frequency for 
    tot_alls = 0  # Count of the total alleles seen
    for cross in sorted(crosses, key=lambda c: c.gen):
        indv_genotypes = site_genotypes[cross.pop]
        for sample in indv_genotypes:
            # Increase sample tally
            info_ns += 1
            genotype = sorted(indv_genotypes[sample])
            # Increase the allele tally
            for al in genotype:
                tot_alls += 1
                if al == 1:
                    info_af += 1
            # Formatted genotype STR
            geno_str = f'{genotype[0]}/{genotype[1]}'
            row.append(geno_str)
    # Add the formatted INFO tab
    info_af = info_af/tot_alls
    info = f'NS={info_ns};AF={info_af:0.04g}'
    row[7] = info
    # Save to the file
    row_str = '\t'.join(row)
    vcf.write(f'{row_str}\n')

def simulate_genotypes(ancestry_map, crosses, parental_allele_freqs, var_site_info, vcf, log=sys.stdout):
    '''Simulate the genotypes for all the possible crosses using discrete parental
    genotyoes. Store the resulting data in the output VCF.'''
    print('\nSimulating genotypes...', file=log, flush=True)
    n_sites = 0
    # First, let's regenerate the VCF header
    make_vcf_header(crosses, ancestry_map, vcf)
    # We process one variant site at a time. They are independent from one another.
    for site in var_site_info:
        n_sites += 1
        # Coordinate and allele info to re-generate the VCF
        info = var_site_info[site]
        # Initialize the genotype dictionary for that given site
        site_genotypes = dict()
        # Now, work through the crosses. They must happen in order of the generations.
        for cross in sorted(crosses, key=lambda c: c.gen):
            # Genotypes for that cross
            site_genotypes.setdefault(cross.pop, dict())
            pop_samples = ancestry_map[cross.pop]
            if not cross.hybrid:
                # Process the parents
                # Get the parental allele freq at that site
                site_freqs = parental_allele_freqs[info.id]
                pop_genotypes = simulate_parental_site_genotypes(pop_samples, cross.pop, site_freqs)
                site_genotypes[cross.pop] = pop_genotypes
            else:
                # Process all the other crossses, based on the fixed, existing genotypes
                pop_genotypes = sample_hybrid_genotypes(pop_samples, cross, site_genotypes)
                site_genotypes[cross.pop] = pop_genotypes
        format_vcf_row(info, site_genotypes, crosses, vcf)
    # Report and Return
    print(f'    Simulated genotypes for {n_sites:,} variant sites.',
          file=log, flush=True)
    vcf.close()


def find_pop_allele_freq(cross, site_freqs):
    '''Find the allele frequencies of a target population based on
    the frequency of the parental population'''
    assert isinstance(cross, Cross)
    p1_freq = site_freqs[cross.par1] # Frequencies of parent1
    p_p1 = p1_freq[0]                # P freq for parent1
    p2_freq = site_freqs[cross.par2] # Frequencies of parent2
    p_p2 = p2_freq[0]                # P freq for parent2
    # Find the new frequency
    p_f1, q_f1 = f1_freq(p_p1, p_p2)
    return [p_f1, q_f1]


def simulate_genotypes_from_prob(ancestry_map, crosses, parental_allele_freqs, var_site_info, vcf, log=sys.stdout):
    '''Simulate the genotypes for all the possible crosses, using the probaility
    (allele frequency) model. Store the resulting data in the output VCF.'''
    print('\nSimulating genotypes...', file=log, flush=True)
    n_sites = 0
    # First, let's regenerate the VCF header
    make_vcf_header(crosses, ancestry_map, vcf)
    # We process one variant site at a time. They are independent from one another.
    for site in var_site_info:
        n_sites += 1
        # Coordinate and allele info to re-generate the VCF
        info = var_site_info[site]
        # Initialize the genotype/allele freq dictionary for that given site
        site_genotypes = dict()
        site_freqs = dict()
        # Now, work through the crosses. They must happen in order of the generations.
        for cross in sorted(crosses, key=lambda c: c.gen):
            # print(cross, site_freqs)
            # Genotypes for that cross
            site_genotypes.setdefault(cross.pop, dict())
            pop_samples = ancestry_map[cross.pop]
            if not cross.hybrid:
                # Process the parents
                # Get the parental allele freq at that site, based on the empirical parental allele frequency
                site_freqs = parental_allele_freqs[info.id]
                # Simulate genotypes from the allele freq
                pop_genotypes = sample_genotypes_from_allele_freqs(pop_samples, cross, site_freqs)
                site_genotypes[cross.pop] = pop_genotypes
            else:
                # Process all the other crossses
                # First, calculate an allele freq for the target population based on the allele freq at the previous generation
                pop_allele_freq = find_pop_allele_freq(cross, site_freqs)
                site_freqs[cross.pop] = pop_allele_freq
                # Sample genotypes based on the calculatated allele frequencies of the previous generation
                pop_genotypes = sample_genotypes_from_allele_freqs(pop_samples, cross, site_freqs)
                site_genotypes[cross.pop] = pop_genotypes
        format_vcf_row(info, site_genotypes, crosses, vcf)
    # Report and Return
    print(f'    Simulated genotypes for {n_sites:,} variant sites.',
          file=log, flush=True)
    vcf.close()

def f1_freq(p_p1, p_p2):
    '''Calculate allele frequencies after a cross following expected Hardy Weinberg proportions'''
    # Check inputs
    assert type(p_p1) is float
    assert 0<=p_p1<=1
    assert type(p_p2) is float
    assert 0<=p_p2<=1
    # q = 1 - p
    q_p1 = 1-p_p1
    q_p2 = 1-p_p2
    # p homozygoys genotype (pp) is p^2
    p1p2 = p_p1*p_p2
    # q homozygoys genotype (qq) is q^2
    q1q2 = q_p1*q_p2
    # pq heterozygoous genotype is 2pq
    pq12 = (p_p1*q_p2)+(p_p2*q_p1)
    # p^2 + 2pq + q^2 = 1
    assert math.isclose((p1p2+q1q2+pq12),1)
    # P(p) = p^2 + pq
    p_f1 = p1p2+(pq12/2)
    # P(q) = q^2 + pq
    q_f1 = q1q2+(pq12/2)
    # p + q = 1
    assert math.isclose((p_f1+q_f1),1)
    return p_f1, q_f1

def main():
    args = parse_args()
    # Initialize log file
    log_f = open(f'{args.outdir}/hybrid_sims.log', 'w')
    initialize_log(args, log_f)
    # Initialize output VCF
    out_vcf = gzip.open(f'{args.outdir}/simulated_hybrids.vcf.gz', 'wt')
    # Get parental IDs and populations
    parent_ids = parse_popmap(args.popmap, log_f)
    # Extract parental alleles from VCF
    parental_allele_freqs, snp_info = parse_vcf(args.vcf, parent_ids, out_vcf, args.min_maf, log_f)
    # Make a map of the simulated crosses and new individuals
    crosses_map, ancestry_map = map_simulated_crosses(args.generations, args.n_individuals, args.outdir, log_f)
    # Simulate the genotypes and save this to a the VCF
    prob = True
    # TODO: For now we are hardcoding to always use allele freqs, but it could be an option in the future. Need to run more tests on this.
    if prob:
        # This function uses a new method that samples genotypes from the allele frequencies
        simulate_genotypes_from_prob(ancestry_map, crosses_map, parental_allele_freqs, snp_info, out_vcf, log_f)
    else:
        # Previous method that sample genotypes from the already existing genotypes in the previous generation
        simulate_genotypes(ancestry_map, crosses_map, parental_allele_freqs, snp_info, out_vcf, log_f) 

    # Close outputs
    log_f.write(f'\n{PROG} finished on {now()}\n')
    log_f.close()

# Run Code
if __name__ == '__main__':
    main()
