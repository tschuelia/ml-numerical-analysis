import sys
sys.path.append(snakemake.scriptdir + "/../../..")

from snakelib.save_trees import save_overall_best_tree

tree_files = snakemake.input.trees
log_files = snakemake.input.logs

out_file = snakemake.output.best_overall_tree

save_overall_best_tree(tree_files, log_files, out_file)