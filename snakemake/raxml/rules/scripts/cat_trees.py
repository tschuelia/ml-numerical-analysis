infiles = snakemake.input.trees
outfile = snakemake.output.best_trees_all_runs

output = []

for tree_file in infiles:
    with open(tree_file) as f:
        content = f.readlines()[0]
        output.append(content.strip())

with open(outfile, "w") as f:
    f.write("\n".join(output))