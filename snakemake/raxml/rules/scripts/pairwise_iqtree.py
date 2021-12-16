import json
import sys
sys.path.append(snakemake.scriptdir + "/../../..")

from snakelib.utils import pairwise_iqtree

msa = snakemake.input.msa
filtered_trees = open(snakemake.input.filtered_trees).readlines()
best_tree = open(snakemake.input.best_tree).readline()
iqtree_command = snakemake.params.iqtree_command
model = snakemake.params.model
max_workers = snakemake.params.max_workers

output = snakemake.output.iqtree_results

all_results = pairwise_iqtree(iqtree_command, filtered_trees, best_tree, msa, model, max_workers)

with open(output, "w") as f:
    json.dump(all_results, f)