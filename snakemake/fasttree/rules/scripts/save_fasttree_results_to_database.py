import sys

sys.path.append(snakemake.scriptdir + "/../../..")

from snakelib import database as db

from fasttree_parser import create_fasttree

db.fasttree_db.init(snakemake.output.database)
db.fasttree_db.connect()
db.fasttree_db.create_tables(
    [
        db.Fasttree,
        db.FasttreeTreesearchTree,
    ]
)

# fmt: off
params_file_paths = snakemake.input.params_file

treesearch_log_file_paths        = snakemake.input.treesearch_log
best_treesearch_tree_file_paths  = snakemake.input.best_treesearch_tree
treesearch_trees_file_paths      = snakemake.input.treesearch_trees
# fmt: on

num_runs = len(params_file_paths)
fasttree_objects = []

for i in range(num_runs):
    # fmt: off
    fasttree = create_fasttree(
        parameter_file_path             = params_file_paths[i],
        treesearch_log_file_path        = treesearch_log_file_paths[i],
        best_tree_file_path             = best_treesearch_tree_file_paths[i],
        all_treesearch_trees_file_path  = treesearch_trees_file_paths[i],
    )

    fasttree_db = db.Fasttree.create(
        blmin   = fasttree.blmin,
        num_pars_trees          = fasttree.num_trees,
        num_rand_trees          = 0,
        best_treesearch_llh     = fasttree.best_treesearch_llh,
        treesearch_total_time   = fasttree.treesearch_total_time,
    )
    # fmt: on

    # FasttreeTreesearchTree
    fasttree.db_best_treesearch_tree_object = None

    for tree_idx in range(fasttree.get_num_of_trees()):
        tree_values = {}
        # fmt: off
        tree_values["llh"]          = fasttree.get_treesearch_llh_for_tree_index(tree_idx)
        tree_values["compute_time"] = fasttree.get_treesearch_compute_time_for_tree_index(tree_idx)
        tree_values["newick_tree"]  = fasttree.get_newick_tree_for_tree_index(tree_idx)

        is_best = (
            fasttree.tree_for_index_is_best(tree_idx)
            and not fasttree.db_best_treesearch_tree_object
        )
        tree_values["is_best"]  = is_best

        tree_values["program"]  = fasttree_db
        tree_values["seed"]     = fasttree.get_treesearch_seed_for_tree_index(tree_idx)
        # fmt: on

        fasttree_treesearch_tree = db.FasttreeTreesearchTree.create(**tree_values)

        if is_best:
            fasttree.db_best_treesearch_tree_object = fasttree_treesearch_tree
