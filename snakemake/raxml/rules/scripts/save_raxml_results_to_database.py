import sys

sys.path.append(snakemake.scriptdir + "/../../..")

import pickle

from snakelib import database as db
from snakelib.database import insert_program_data, insert_treesarch_data, insert_eval_data, insert_statstest_data
from snakelib.utils import get_iqtree_results_for_eval_tree_str
from snakelib.iqtree_statstest_parser import get_iqtree_results

from raxml_parser import create_raxml

db.raxml_db.init(snakemake.output.database)
db.raxml_db.connect()
db.raxml_db.create_tables(
    [
        db.Raxmlng,
        db.RaxmlTreesearchTree,
        db.RaxmlEvalTree,
        db.RaxmlEvalTreeStatsTest,
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

filtered_trees_clusters     = snakemake.input.filtered_trees_clusters
iqtree_statstests_results   = get_iqtree_results(snakemake.input.iqtree_statstests_results)

# fmt: on
num_runs = len(params_file_paths)
eval_tree_objects = []

for i in range(num_runs):
    # fmt:off
    raxml = create_raxml(
        parameter_file_path             = params_file_paths[i],
        treesearch_log_file_path        = treesearch_log_file_paths[i],
        eval_log_file_path              = eval_log_file_paths[i],
        best_tree_file_path             = best_treesearch_tree_file_paths[i],
        all_treesearch_trees_file_path  = treesearch_trees_file_paths[i],
        best_eval_tree_file_path        = best_eval_tree_log_file_paths[i],
        raxml_command                   = snakemake.params.raxml_command,
        all_eval_trees_file_path        = eval_trees_file_paths[i],
    )
    # fmt: on
    raxml_db = insert_program_data(raxml, db.Raxmlng)

    # RaxmlTreesearchTree
    _, best_tree = insert_treesarch_data(raxml, raxml_db, db.RaxmlTreesearchTree)
    raxml.db_best_treesearch_tree_object = best_tree

    # RaxmlEvalTree for best RaxmlTreesearchTree (raxml.db_best_treesearch_tree_object)
    eval_trees, best_eval_tree = insert_eval_data(raxml, raxml.db_best_treesearch_tree_object, db.RaxmlEvalTree)
    eval_tree_objects += eval_trees
    raxml.db_best_eval_tree = best_eval_tree

# Iqtree significance tests
with open(filtered_trees_clusters, "rb") as f:
    clusters = pickle.load(f)

all_results = []
cluster_ids = []

for tree in eval_tree_objects:
    results, cluster_id = get_iqtree_results_for_eval_tree_str(iqtree_statstests_results, tree.newick_tree, clusters)
    all_results.append(results)
    cluster_ids.append(cluster_id)

insert_statstest_data(
    eval_trees=eval_tree_objects,
    statstest_results=all_results,
    cluster_ids=cluster_ids,
    database_table=db.RaxmlEvalTreeStatsTest
)
