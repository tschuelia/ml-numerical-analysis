import sys

sys.path.append(snakemake.scriptdir + "/../../..")
from snakelib.save_best_tree import save_best_tree

tree_file = snakemake.input.trees
logs_file = snakemake.input.logs
output_file = snakemake.output.best_tree
program_name = snakemake.params.program

save_best_tree(tree_file, logs_file, output_file, program_name)