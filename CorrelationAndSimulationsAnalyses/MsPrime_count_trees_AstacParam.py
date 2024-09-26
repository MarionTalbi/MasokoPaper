################################################################
####   Count the number of recombination event after split  ####
################################################################

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

###Arguments
argparser.add_argument("-r", "--recombination_map", type = str, help="Path to a recombination map file") 
argparser.add_argument("-t", "--split_time", type=float, help="Split time between the populations (in generations)",default=2484)
argparser.add_argument("-u", "--mut_rate", type=float, help="Mutation rate per bp per generation",default=3.5e-9)
argparser.add_argument("-n", "--num_indiv", help="The number of (diploid) individuals to simulate from each population/species", type=int, default=35)
argparser.add_argument("-s", "--seed", help="The random seed to run the msprime simulation", type=int, default=123456)
#argparser.add_argument("-f", "--fullARG", help="Record the full ARG in msprime simulation (still experimental)", type=bool, default=False)
args = argparser.parse_args()


##Add here : to be able to implement a non homogeneous recombination map:
## Read in a recombination map A:
with open(args.recombination_map, "r") as cichlid_map_file:
   rate_map = msprime.RateMap.read_hapmap(cichlid_map_file)

split_time = args.split_time

##Implementation of a more complex demographic history :
demography = msprime.Demography()
demography.add_population(name="Ben", initial_size=72630)
demography.add_population(name="Lit", initial_size=75996)
demography.add_population(name="Common", initial_size=12907)
demography.set_migration_rate(source="Ben", dest="Lit", rate=0.000115) ##Low migration rates
demography.set_migration_rate(source="Lit", dest="Ben", rate=0.0000701) ##Low migration rates
demography.add_population_parameters_change(time=1943, population="Lit", initial_size=25195)
demography.add_population_parameters_change(time=2006, population="Ben", initial_size=10295)
demography.add_population_split(time=2484, derived=["Ben", "Lit"], ancestral="Common")
demography.add_population_parameters_change(time=2880, population="Common", initial_size=5699)
demography.add_population_parameters_change(time=4134, population="Common", initial_size=6792)
demography.add_population_parameters_change(time=7441, population="Common", initial_size=21764)
demography.add_population_parameters_change(time=18379, population="Common", initial_size=84371)
demography.add_population_parameters_change(time=54395, population="Common", initial_size=177574)
demography.add_population_parameters_change(time=112126, population="Common", initial_size=103772)
demography.debug()

seed = args.seed
tsA = msprime.sim_ancestry(
    {"Ben": args.num_indiv},
    demography=demography,
    recombination_rate=rate_map,
    random_seed=seed,
    end_time=split_time,
)

tsB = msprime.sim_ancestry(
    {"Lit": args.num_indiv},
    demography=demography,
    recombination_rate=rate_map,
    random_seed=seed + 1,
    end_time=split_time,
)

tablesA = tsA.dump_tables() # Returns a modifiable copy of the tables defining this tree sequence.
tablesB = tsB.dump_tables() # SAme (table are : individuals, nodes, edges, ..)

tablesB.edges.child += tsA.num_nodes ## tablesB.edges.child = tablesB.edges.child + tsA.num_nodes
tablesB.edges.parent += tsA.num_nodes
tablesB.nodes.individual += tsA.num_individuals # Add this time the num_individuals

individuals_dict = tablesB.individuals.asdict()
del individuals_dict["metadata_schema"]
tablesA.individuals.append_columns(**individuals_dict)
nodes_dict = tablesB.nodes.asdict()
del nodes_dict["metadata_schema"]
tablesA.nodes.append_columns(**nodes_dict)
edges_dict = tablesB.edges.asdict()
del edges_dict["metadata_schema"]
tablesA.edges.append_columns(**edges_dict)
tablesA.sort()
tablesA.build_index()
ts_joint = tablesA.tree_sequence()

tsC = msprime.sim_ancestry(
    initial_state=ts_joint,
    demography=demography,
    recombination_rate=rate_map,
    random_seed=seed + 2,
)


print('Proportion of recombination events after split (based on "num_trees"):')
print( ( (tsA.num_trees-1) + (tsB.num_trees-1) ) / (tsC.num_trees-1) )
