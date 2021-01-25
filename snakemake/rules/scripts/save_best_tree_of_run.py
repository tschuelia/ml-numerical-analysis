from utils import get_all_raxml_llhs

with open(snakemake.input.all_trees) as f:
    all_trees = f.readlines()

all_llhs = get_all_raxml_llhs(snakemake.input.all_logs)
idx_best = all_llhs.index(max(all_llhs))

with open(snakemake.output.best_tree_of_run, "w") as f:
    f.write(all_trees[idx_best])