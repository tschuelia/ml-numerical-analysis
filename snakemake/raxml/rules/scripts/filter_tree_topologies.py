import sys
sys.path.append(snakemake.scriptdir + "/../../..")

from snakelib.utils import filter_tree_topologies, read_file_contents

eval_trees = read_file_contents(snakemake.input.all_trees)

filter_tree_topologies(
    raxml_command=snakemake.params.raxml_command,
    eval_trees=eval_trees,
    log_path=snakemake.output.raxml_rfdist_log,
    filtered_trees_path=snakemake.output.filtered_trees,
    clusters_path=snakemake.output.clusters,
)