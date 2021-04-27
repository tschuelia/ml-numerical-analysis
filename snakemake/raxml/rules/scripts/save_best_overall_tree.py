import json
import pandas as pd

tree_files = snakemake.input.trees
log_files = snakemake.input.logs

out_file = snakemake.output.best_overall_tree

trees = []
logs = []

for tf in tree_files:
    with open(tf) as f:
        trees.append(f.readlines()[0])

for l in log_files:
    with open(l) as f:
        logs.append(json.load(f))

df = pd.DataFrame(logs)
df["newick"] = trees
df = df.sort_values(by=["llh", "runtime"], ascending=[False, True])

with open(out_file, "w") as f:
    f.write(df.head(1)["newick"].item())