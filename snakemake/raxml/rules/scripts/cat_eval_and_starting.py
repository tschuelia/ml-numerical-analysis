#"cat {input.starting_eval_trees} {input.eval_trees} > {output.all_trees}"

starting_eval = snakemake.input.starting_eval_trees
eval = snakemake.input.eval_trees

output = open(starting_eval).read() + "\n" + open(eval).read()

with open(snakemake.output.all_trees, "w") as f:
    f.write(output)
