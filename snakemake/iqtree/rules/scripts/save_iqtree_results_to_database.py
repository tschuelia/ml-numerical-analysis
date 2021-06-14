import sys

sys.path.append(snakemake.scriptdir + "/../../..")

from peewee import chunked

from snakelib import database as db

from iqtree_parser import create_iqtree, create_Experiment

db.iqtree_db.init(snakemake.output.database)
db.iqtree_db.connect()
db.iqtree_db.create_tables(
    [
        db.Iqtree,
        db.IqtreeTreesearchTree,
        db.IqtreeEvalTree,
        db.IqtreeEvalTreeStatsTest
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

# all_best_eval_trees_file_paths            = snakemake.input.all_best_eval_trees
# best_overall_eval_tree_file_paths      = snakemake.input.best_overall_eval_tree
# iqtree_significance_summary_file_paths = snakemake.input.iqtree_significance_summary
# fmt: on

num_runs = len(params_file_paths)
iqtree_objects = []
best_eval_tree_objects = []

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
        lh_eps  = iqtree.model_param_epsilon,
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
        is_best = iqtree.eval_tree_for_index_is_best(eval_tree_idx)
        # fmt: off
        eval_tree_values["start_tree"]  = iqtree.db_best_treesearch_tree_object
        eval_tree_values["llh"]         = iqtree.get_eval_llh_for_tree_index(eval_tree_idx)
        eval_tree_values["newick_tree"] = iqtree.get_newick_eval_tree_for_tree_index(eval_tree_idx)

        eval_tree_values["compute_time"] = iqtree.get_eval_compute_time_for_tree_index(eval_tree_idx)

        eval_tree_values["is_best"]     = is_best

        eval_tree_values["number_of_taxa"]          = iqtree.get_number_of_taxa_for_eval_tree_index(eval_tree_idx)
        eval_tree_values["total_branch_length"]     = iqtree.get_total_branch_length_for_eval_tree_index(eval_tree_idx)
        eval_tree_values["average_branch_length"]   = iqtree.get_average_branch_length_for_eval_tree_index(eval_tree_idx)

        eval_tree_values["eval_blmin"]  = iqtree.get_eval_blmin_for_tree_index(eval_tree_idx)
        eval_tree_values["eval_blmax"]  = iqtree.get_eval_blmax_for_tree_index(eval_tree_idx)
        eval_tree_values["eval_lh_eps"] = iqtree.get_eval_model_param_epsilon_for_tree_index(eval_tree_idx)
        # fmt: on
        iqtree_eval_tree = db.IqtreeEvalTree.create(**eval_tree_values)

        if is_best:
            best_eval_tree_objects.append(iqtree_eval_tree)

# # Iqtree significance tests
# # fmt: off
# experiment = create_Experiment(
#     best_trees_file_path                = all_best_eval_trees_file_paths,
#     best_overall_eval_tree_file_path    = best_overall_eval_tree_file_paths,
#     iqtree_statstest_results_file_path  = iqtree_significance_summary_file_paths
# )
# # fmt: on
#
# best_overall_eval_tree = [tree for tree in best_eval_tree_objects if experiment.eval_tree_is_overall_best(tree.newick_tree)]
# assert len(best_overall_eval_tree) >= 1, "Overall best eval tree not in all eval trees."
# best_overall_eval_tree = best_overall_eval_tree[0]
#
# insert_into_significance_table = []
#
# for eval_tree in best_eval_tree_objects:
#     newick_str = eval_tree.newick_tree
#     statstest_values = {}
#     statstest_values["reference_tree_id"] = best_overall_eval_tree
#     statstest_values["tree_id"] = eval_tree
#     statstest_values["iqtree_llh"] = experiment.get_iqtree_llh_for_eval_tree(newick_str)
#     statstest_values["deltaL"] = experiment.get_iqtree_deltaL_for_eval_tree(newick_str)
#
#     tests = experiment.get_iqtree_test_results_for_eval_tree(newick_str)
#
#     if "bp-RELL" in tests:
#         statstest_values["bpRell"] = tests["bp-RELL"]["score"]
#         statstest_values["bpRell_significant"] = tests["bp-RELL"]["significant"]
#
#     if "p-KH" in tests:
#         statstest_values["pKH"] = tests["p-KH"]["score"]
#         statstest_values["pKH_significant"] = tests["p-KH"]["significant"]
#
#     if "p-SH" in tests:
#         statstest_values["pSH"] = tests["p-SH"]["score"]
#         statstest_values["pSH_significant"] = tests["p-SH"]["significant"]
#
#     if "p-WKH" in tests:
#         statstest_values["pWKH"] = tests["p-WKH"]["score"]
#         statstest_values["pWKH_significant"] = tests["p-WKH"]["significant"]
#
#     if "p-WSH" in tests:
#         statstest_values["pWSH"] = tests["p-WSH"]["score"]
#         statstest_values["pWSH_significant"] = tests["p-WSH"]["significant"]
#
#     if "c-ELW" in tests:
#         statstest_values["cELW"] = tests["c-ELW"]["score"]
#         statstest_values["cELW_significant"] = tests["c-ELW"]["significant"]
#
#     if "p-AU" in tests:
#         statstest_values["pAU"] = tests["p-AU"]["score"]
#         statstest_values["pAU_significant"] = tests["p-AU"]["significant"]
#
#     insert_into_significance_table.append(statstest_values)
#
# with db.iqtree_db.atomic():
#     for batch in chunked(insert_into_significance_table, 100):
#         db.IqtreeEvalTreeStatsTest.insert_many(batch).execute()
