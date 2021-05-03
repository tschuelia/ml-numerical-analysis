import sys
sys.path.append(snakemake.scriptdir + "/../../..")

from snakelib.utils import cat_input_to_output

cat_input_to_output(
    input_files=snakemake.input.trees,
    output_file=snakemake.output.best_trees_all_runs
)