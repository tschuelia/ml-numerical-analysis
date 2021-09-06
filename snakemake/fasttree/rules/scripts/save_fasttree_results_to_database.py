import sys

sys.path.append(snakemake.scriptdir + "/../../..")

from peewee import chunked

from snakelib import database as db
from snakelib.database import insert_program_data, insert_treesarch_data

from fasttree_parser import create_fasttree, create_Experiment

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

# best_trees_file_path                   = snakemake.input.best_trees
# best_overall_eval_tree_file_paths      = snakemake.input.best_overall_eval_tree
# iqtree_significance_summary_file_paths = snakemake.input.iqtree_significance_summary

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
    best_tree_objects += [t for t in tree_objects if t.is_best]

# # Iqtree significance tests
# # fmt: off
# experiment = create_Experiment(
#     best_trees_file_path                = best_trees_file_path,
#     best_overall_eval_tree_file_path    = best_overall_eval_tree_file_paths,
#     iqtree_statstest_results_file_path  = iqtree_significance_summary_file_paths
# )
# # fmt: on
#
# best_overall_eval_tree = [tree for tree in best_tree_objects if experiment.eval_tree_is_overall_best(tree.newick_tree)]
# assert len(best_overall_eval_tree) >= 1, "Overall best eval tree not in all eval trees."
# best_overall_eval_tree = best_overall_eval_tree[0]
#
# insert_into_significance_table = []
#
# for tree in best_tree_objects:
#     newick_str = tree.newick_tree
#     statstest_values = {}
#     statstest_values["reference_tree_id"] = best_overall_eval_tree
#     statstest_values["tree_id"] = tree
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
# with db.fasttree_db.atomic():
#     for batch in chunked(insert_into_significance_table, 100):
#         db.FasttreeEvalTreeStatsTest.insert_many(batch).execute()
