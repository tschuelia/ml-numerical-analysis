import sys

sys.path.append(snakemake.scriptdir + "/../../..")

from snakelib import database as db

from iqtree_parser import create_iqtree

db.iqtree_db.init(snakemake.output.database)
db.iqtree_db.connect()
db.iqtree_db.create_tables(
    [
        db.Iqtree,
        db.IqtreeTreesearchTree,
        db.IqtreeEvalTree,
    ]
)

# fmt: off
params_file_paths = snakemake.input.params_file

treesearch_log_file_paths        = snakemake.input.treesearch_log
best_treesearch_tree_file_paths  = snakemake.input.best_treesearch_tree
treesearch_trees_file_paths      = snakemake.input.treesearch_trees

eval_log_file_paths          = snakemake.input.eval_log
best_eval_tree_file_paths    = snakemake.input.best_eval_tree
eval_trees_file_paths        = snakemake.input.eval_trees
# fmt: on

num_runs = len(params_file_paths)
iqtree_objects = []

for i in range(num_runs):
    # fmt: off
    iqtree = create_iqtree(
        parameter_file_path             = params_file_paths[i],
        treesearch_log_file_path        = treesearch_log_file_paths[i],
        eval_log_file_path              = eval_log_file_paths[i],
        best_tree_file_path             = best_treesearch_tree_file_paths[i],
        all_treesearch_trees_file_path  = treesearch_trees_file_paths[i],
        best_eval_tree_file_path        = best_eval_tree_file_paths[i],
        all_eval_trees_file_path        = eval_trees_file_paths[i],
    )

    iqtree_db = db.Iqtree.create(
        blmin   = iqtree.blmin,
        blmax   = iqtree.blmax,
        lh_eps  = iqtree.lh_eps,
        num_pars_trees          = iqtree.num_pars_trees,
        num_rand_trees          = 0, #iqtree.num_rand_trees,
        best_treesearch_llh     = iqtree.best_treesearch_llh,
        best_evaluation_llh     = iqtree.best_evaluation_llh,
        treesearch_total_time   = iqtree.treesearch_total_time,
    )
    # fmt: on

    # IqtreeTreesearchTree
    iqtree.db_best_treesearch_tree_object = None

    for tree_idx in range(iqtree.get_num_of_trees()):
        tree_values = {}
        # fmt: off
        tree_values["llh"]          = iqtree.get_treesearch_llh_for_tree_index(tree_idx)
        tree_values["compute_time"] = iqtree.get_treesearch_compute_time_for_tree_index(tree_idx)
        tree_values["newick_tree"]  = iqtree.get_newick_tree_for_tree_index(tree_idx)

        is_best = (
            iqtree.tree_for_index_is_best(tree_idx)
            and not iqtree.db_best_treesearch_tree_object
        )
        tree_values["is_best"]  = is_best
        tree_values["number_of_taxa"]       = iqtree.get_number_of_taxa_for_tree_index(tree_idx)
        tree_values["total_branch_length"]  = iqtree.get_total_branch_length_for_tree_index(tree_idx)
        tree_values["average_branch_length"] = iqtree.get_average_branch_length_for_tree_index(tree_idx)

        tree_values["program"]  = iqtree_db
        tree_values["seed"]     = iqtree.get_treesearch_seed_for_tree_index(tree_idx)
        # fmt: on

        iqtree_treesearch_tree = db.IqtreeTreesearchTree.create(**tree_values)

        if is_best:
            iqtree.db_best_treesearch_tree_object = iqtree_treesearch_tree

    # IqtreeEvalTree for best IqtreeTreesearchTree (iqtree.db_best_treesearch_tree_object)
    for eval_tree_idx in range(iqtree.get_num_of_eval_trees()):
        eval_tree_values = {}
        # fmt: off
        eval_tree_values["start_tree"]  = iqtree.db_best_treesearch_tree_object
        eval_tree_values["llh"]         = iqtree.get_eval_llh_for_tree_index(eval_tree_idx)
        eval_tree_values["newick_tree"] = iqtree.get_newick_eval_tree_for_tree_index(eval_tree_idx)

        eval_tree_values["compute_time"] = iqtree.get_eval_compute_time_for_tree_index(eval_tree_idx)

        eval_tree_values["is_best"]     = iqtree.eval_tree_for_index_is_best(eval_tree_idx)

        eval_tree_values["number_of_taxa"]          = iqtree.get_number_of_taxa_for_eval_tree_index(eval_tree_idx)
        eval_tree_values["total_branch_length"]     = iqtree.get_total_branch_length_for_eval_tree_index(eval_tree_idx)
        eval_tree_values["average_branch_length"]   = iqtree.get_average_branch_length_for_eval_tree_index(eval_tree_idx)

        eval_tree_values["eval_blmin"]  = iqtree.get_eval_blmin_for_tree_index(eval_tree_idx)
        eval_tree_values["eval_blmax"]  = iqtree.get_eval_blmax_for_tree_index(eval_tree_idx)
        eval_tree_values["eval_lh_eps"] = iqtree.get_eval_lh_eps_for_tree_index(eval_tree_idx)
        # fmt: on
        iqtree_eval_tree = db.IqtreeEvalTree.create(**eval_tree_values)
