############################################################
######              MSPRIME SIMULATIONS               ######
############################################################
import sys
if sys.version_info[0] < 3:
    print('Python 3 is needed to run this script!')
    sys.exit(0)
import argparse
import msprime
import tskit
import re

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    description="Simulate SNPs in a VCF format using msprime with a custom recombination map. ", add_help=True)

argparser.add_argument("-r", "--recombination_map", type = str, help="Path to a recombination map file") 
argparser.add_argument("-c", "--chromosome", type = str, help="Name of the chr")
argparser.add_argument("-t", "--split_time", type=float, help="Split time between the populations (in generations)",default=2484)
argparser.add_argument("-u", "--mut_rate", type=float, help="Mutation rate per bp per generation",default=3.5e-9)
argparser.add_argument("-n", "--num_indiv", help="The number of (diploid) individuals to simulate from each population/species", type=int, default=35)
argparser.add_argument("-s", "--seed", help="The random seed to run the msprime simulation", type=int, default=123456)
args = argparser.parse_args()

with open(args.recombination_map, "r") as cichlid_map_file:
   rate_map = msprime.RateMap.read_hapmap(cichlid_map_file)

split_time = args.split_time

demography = msprime.Demography()
demography.add_population(name="Ben", initial_size=72630)
demography.add_population(name="Lit", initial_size=75996)
demography.add_population(name="Common", initial_size=12907)
demography.set_migration_rate(source="Ben", dest="Lit", rate=0.000115) ##Migration rate can be modify here
demography.set_migration_rate(source="Lit", dest="Ben", rate=0.0000701) ##Migration rate can be modify here
demography.add_population_parameters_change(time=1943, population="Lit", initial_size=25195)
demography.add_population_parameters_change(time=2006, population="Ben", initial_size=10295)
demography.add_population_split(time=split_time, derived=["Ben", "Lit"], ancestral="Common")
demography.add_population_parameters_change(time=2880, population="Common", initial_size=5699)
demography.add_population_parameters_change(time=4134, population="Common", initial_size=6792)
demography.add_population_parameters_change(time=7441, population="Common", initial_size=21764)
demography.add_population_parameters_change(time=18379, population="Common", initial_size=84371)
demography.add_population_parameters_change(time=54395, population="Common", initial_size=177574)
demography.add_population_parameters_change(time=112126, population="Common", initial_size=103772)
demography.debug()

ts = msprime.sim_ancestry(
    {"Lit": args.num_indiv,"Ben": args.num_indiv},
    demography=demography,
    recombination_rate=rate_map,
    random_seed=args.seed
)

mutated_ts = msprime.sim_mutations(ts, rate=args.mut_rate)

with open("./vcf/Msprime_" + str(args.chromosome) + "_simul.vcf", "w") as vcf_file:
        mutated_ts.write_vcf(vcf_file)
