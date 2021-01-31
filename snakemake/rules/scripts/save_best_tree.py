from utils import get_all_raxml_llhs, get_all_iqtree_llhs

with open(snakemake.input.trees) as f:
    all_trees = f.readlines()

if snakemake.params.program == "raxml":
    all_llhs = get_all_raxml_llhs(snakemake.input.logs)
elif snakemake.params.program == "iqtree":
    all_llhs = get_all_iqtree_llhs(snakemake.input.logs)
else:
    raise ValueError(
        "Either the specified tree inference program is not supported or you did not specify any."
    )

idx_best = all_llhs.index(max(all_llhs))

with open(snakemake.output.best_tree, "w") as f:
    f.write(all_trees[idx_best])