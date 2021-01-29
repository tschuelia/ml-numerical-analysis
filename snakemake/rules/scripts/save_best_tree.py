from utils import get_all_raxml_llhs

with open(snakemake.input.trees) as f:
    all_trees = f.readlines()

all_llhs = get_all_raxml_llhs(snakemake.input.logs)
idx_best = all_llhs.index(max(all_llhs))

with open(snakemake.output.best_tree, "w") as f:
    f.write(all_trees[idx_best])