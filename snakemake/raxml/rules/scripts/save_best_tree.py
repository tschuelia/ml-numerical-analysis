import sys
sys.path.append(snakemake.scriptdir + "/../../..")

from raxml_utils import get_all_raxml_llhs, get_raxml_elapsed_time
from snakelib.save_trees import save_best_tree_and_log

tree_file = snakemake.input.trees
logs_file = snakemake.input.logs
tree_output_file = snakemake.output.best_tree
log_output_file = snakemake.output.best_log

all_llhs = get_all_raxml_llhs(logs_file)
all_runtimes = get_raxml_elapsed_time(logs_file)

save_best_tree_and_log(
    tree_file=tree_file,
    all_llhs=all_llhs,
    all_runtimes=all_runtimes,
    out_tree=tree_output_file,
    out_log=log_output_file
)