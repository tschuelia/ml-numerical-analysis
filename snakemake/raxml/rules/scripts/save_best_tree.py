import sys
import json

sys.path.append(snakemake.scriptdir + "/../../..")

from raxml_utils import get_all_raxml_llhs, get_raxml_elapsed_time

tree_file = snakemake.input.trees
logs_file = snakemake.input.logs
tree_output_file = snakemake.output.best_tree
log_output_file = snakemake.output.best_log

with open(tree_file) as f:
    all_trees = f.readlines()

all_llhs = get_all_raxml_llhs(logs_file)
all_runtimes = get_raxml_elapsed_time(logs_file)

idx_best = all_llhs.index(max(all_llhs))

with open(tree_output_file, "w") as f:
    f.write(all_trees[idx_best])

data = {
    "llh": all_llhs[idx_best],
    "runtime": all_runtimes[idx_best]
}

with open(log_output_file, "w") as f:
    json.dump(data, f)