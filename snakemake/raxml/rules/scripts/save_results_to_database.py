import sys

sys.path.append(snakemake.scriptdir + "/../../..")

from snakelib import database as db
from .raxml_parser import create_raxml


db.db.init(snakemake.output.database)
db.db.connect()
db.db.create_tables(
    [
        db.Raxmlng,
        db.RaxmlTreesearchTree,
        db.RaxmlEvalTree,
        db.RFDistTreesearchTree,
        db.RFDistEvalTree,
    ]
)

# fmt: off
params_file_paths = snakemake.input.params_file

treesearch_log_file_paths         = snakemake.input.treesearch_log
best_treesearch_tree_file_paths   = snakemake.input.best_treesearch_tree
treesearch_trees_file_paths       = snakemake.input.treesearch_trees

eval_log_file_paths           = snakemake.input.eval_log
best_eval_tree_log_file_paths = snakemake.input.best_eval_tree
eval_trees_file_paths         = snakemake.input.eval_trees

#iqtree_statstest_results_file_paths = snakemake.input.iqtree_statstest_results

treesearch_rfDist_log_file_paths  = snakemake.input.treesearch_rfDist_log
treesearch_rfDist_file_paths      = snakemake.input.treesearch_rfDist

all_best_treesearch_trees_file_paths      = snakemake.input.all_best_treesearch_trees[0]
rfDist_best_treesearch_trees_file_paths   = snakemake.input.rfDist_best_treesearch_trees[0]

all_best_eval_trees_file_paths            = snakemake.input.all_best_eval_trees[0]
rfDist_best_eval_trees_file_paths         = snakemake.input.rfDist_best_eval_trees[0]

# fmt: on

num_runs = len(params_file_paths)