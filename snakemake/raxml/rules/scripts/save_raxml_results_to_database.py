import sys

sys.path.append(snakemake.scriptdir + "/../../..")

from peewee import chunked

from snakelib import database as db

from raxml_parser import create_raxml, create_Experiment

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

#best_overall_eval_tree_file_paths      = snakemake.input.best_overall_eval_tree
#iqtree_significance_summary_file_paths = snakemake.input.iqtree_significance_summary

# fmt: on

num_runs = len(params_file_paths)
raxml_objects = []
best_eval_tree_objects = []

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
    raxml_objects.append(raxml)

    # fmt:off
    raxml_db = db.Raxmlng.create(
        blmin   = raxml.blmin,
        blmax   = raxml.blmax,
        lh_eps  = raxml.lh_epsilon,
        raxml_param_epsilon     = raxml.model_param_epsilon,
        branch_length_smoothing = raxml.branch_length_smoothing,
        spr_lh_epsilon          = raxml.spr_lh_epsilon,
        bfgs_factor             = raxml.bfgs_factor,

        num_pars_trees          = raxml.num_pars_trees,
        num_rand_trees          = raxml.num_rand_trees,
        best_treesearch_llh     = raxml.best_treesearch_llh,
        best_evaluation_llh     = raxml.best_evaluation_llh,
        treesearch_total_time   = raxml.treesearch_total_time,
    )
    # fmt: on

    # RaxmlTreesearchTree
    raxml_treesearch_tree_db_objects = []
    raxml.db_best_treesearch_tree_object = None

    for tree_idx in range(raxml.get_number_of_trees()):
        tree_values = {}
        tree_values["llh"] = raxml.treesearch_llhs[tree_idx]
        tree_values["compute_time"] = raxml.treesearch_compute_times[tree_idx]
        tree_values["newick_tree"] = raxml.treesearch_trees[tree_idx].newick_str

        is_best = (
            raxml.tree_for_index_is_best(tree_idx)
            and not raxml.db_best_treesearch_tree_object
        )
        tree_values["is_best"] = is_best
        tree_values["number_of_taxa"] = raxml.treesearch_trees[tree_idx].number_of_taxa
        tree_values["total_branch_length"] = raxml.treesearch_trees[tree_idx].total_branch_length
        tree_values["average_branch_length"] = raxml.treesearch_trees[tree_idx].average_branch_length

        tree_values["program"] = raxml_db
        tree_values["seed"] = raxml.treeseach_seeds[tree_idx]

        raxml_treesearch_tree = db.RaxmlTreesearchTree.create(**tree_values)
        raxml_treesearch_tree_db_objects.append(raxml_treesearch_tree)

        if is_best:
            raxml.db_best_treesearch_tree_object = raxml_treesearch_tree


    # RaxmlEvalTree for best RaxmlTreesearchTree (raxml.db_best_treesearch_tree_object)
    for eval_tree_idx in range(raxml.get_number_of_eval_trees()):
        is_best = raxml.eval_tree_for_index_is_best(eval_tree_idx)
        # fmt: off
        eval_tree_values = {}
        eval_tree_values["start_tree"]                  = raxml.db_best_treesearch_tree_object
        eval_tree_values["llh"]                         = raxml.eval_llhs[eval_tree_idx]
        eval_tree_values["newick_tree"]                 = raxml.eval_trees[eval_tree_idx].newick_str
        eval_tree_values["compute_time"]                = raxml.eval_compute_times[eval_tree_idx]
        eval_tree_values["is_best"]                     = is_best
        eval_tree_values["number_of_taxa"]              = raxml.eval_trees[eval_tree_idx].number_of_taxa
        eval_tree_values["total_branch_length"]         = raxml.eval_trees[eval_tree_idx].total_branch_length
        eval_tree_values["average_branch_length"]       = raxml.eval_trees[eval_tree_idx].average_branch_length
        eval_tree_values["eval_blmin"]                  = raxml.eval_blmins[eval_tree_idx]
        eval_tree_values["eval_blmax"]                  = raxml.eval_blmaxs[eval_tree_idx]
        eval_tree_values["eval_lh_eps"]                 = raxml.eval_lh_epsilons[eval_tree_idx]
        eval_tree_values["eval_raxml_param_epsilon"]    = raxml.eval_model_param_epsilons[eval_tree_idx]
        eval_tree_values["eval_raxml_brlen_smoothings"] = raxml.eval_raxml_brlen_smoothings[eval_tree_idx]
        eval_tree_values["eval_spr_lh_epsilon"]         = raxml.eval_spr_lh_epsilons[eval_tree_idx]
        eval_tree_values["eval_bfgs_factor"]            = raxml.eval_bfgs_factors[eval_tree_idx]

        # fmt: on
        raxml_eval_tree = db.RaxmlEvalTree.create(**eval_tree_values)

        if is_best:
            raxml.db_best_eval_tree = raxml_eval_tree
            best_eval_tree_objects.append(raxml_eval_tree)


# Iqtree significance tests
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
# with db.raxml_db.atomic():
#     for batch in chunked(insert_into_significance_table, 100):
#        db.RaxmlEvalTreeStatsTest.insert_many(batch).execute()
