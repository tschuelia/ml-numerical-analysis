import sys

sys.path.append(snakemake.scriptdir + "/../../..")

import pickle

from snakelib import database as db
from snakelib.database import insert_program_data, insert_treesarch_data, insert_statstest_data
from snakelib.utils import get_iqtree_results_for_eval_tree_str
from snakelib.iqtree_statstest_parser import get_iqtree_results

from fasttree_parser import create_fasttree

db.fasttree_db.init(snakemake.output.database)
db.fasttree_db.connect()
db.fasttree_db.create_tables(
    [
        db.Fasttree,
        db.FasttreeTreesearchTree,
        db.FasttreeEvalTreeStatsTest,
    ]
)

# fmt: off
params_file_paths = snakemake.input.params_file

treesearch_log_file_paths        = snakemake.input.treesearch_log
best_treesearch_tree_file_paths  = snakemake.input.best_treesearch_tree
treesearch_trees_file_paths      = snakemake.input.treesearch_trees

filtered_trees_clusters     = snakemake.input.filtered_trees_clusters
iqtree_statstests_results   = get_iqtree_results(snakemake.input.iqtree_statstests_results)
# fmt: on

num_runs = len(params_file_paths)
fasttree_objects = []
best_tree_objects = []

for i in range(num_runs):
    # fmt: off
    fasttree = create_fasttree(
        parameter_file_path             = params_file_paths[i],
        treesearch_log_file_path        = treesearch_log_file_paths[i],
        best_tree_file_path             = best_treesearch_tree_file_paths[i],
        all_treesearch_trees_file_path  = treesearch_trees_file_paths[i],
    )

    # fmt: on

    fasttree_db = insert_program_data(fasttree, db.Fasttree)

    # FasttreeTreesearchTree
    tree_objects, best_tree = insert_treesarch_data(fasttree, fasttree_db, db.FasttreeTreesearchTree)
    fasttree.db_best_treesearch_tree_object = best_tree
    best_tree_objects += tree_objects #[t for t in tree_objects if t.is_best]

with open(filtered_trees_clusters, "rb") as f:
    clusters = pickle.load(f)

all_results = []
cluster_ids = []

for tree in best_tree_objects:
    results, cluster_id = get_iqtree_results_for_eval_tree_str(iqtree_statstests_results, tree.newick_tree, clusters)
    all_results.append(results)
    cluster_ids.append(cluster_id)

insert_statstest_data(
    eval_trees=best_tree_objects,
    statstest_results=all_results,
    cluster_ids=cluster_ids,
    database_table=db.FasttreeEvalTreeStatsTest
)