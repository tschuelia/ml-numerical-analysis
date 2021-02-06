import sys

sys.path.append(snakemake.scriptdir + "/../../..")

from raxml_utils import get_all_raxml_llhs

tree_file = snakemake.input.trees
logs_file = snakemake.input.logs
output_file = snakemake.output.best_tree

with open(tree_file) as f:
    all_trees = f.readlines()

all_llhs = get_all_raxml_llhs(logs_file)

idx_best = all_llhs.index(max(all_llhs))

with open(output_file, "w") as f:
    f.write(all_trees[idx_best])
