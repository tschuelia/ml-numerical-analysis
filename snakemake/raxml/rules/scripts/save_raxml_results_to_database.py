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
        db.RFDistTreesearchTree,
        db.RFDistEvalTree,
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

treesearch_rfDist_log_file_paths  = snakemake.input.treesearch_rfDist_log
treesearch_rfDist_file_paths      = snakemake.input.treesearch_rfDist

all_best_treesearch_trees_file_paths      = snakemake.input.all_best_treesearch_trees
rfDist_best_treesearch_trees_file_paths   = snakemake.input.rfDist_best_treesearch_trees

all_best_eval_trees_file_paths            = snakemake.input.all_best_eval_trees
rfDist_best_eval_trees_file_paths         = snakemake.input.rfDist_best_eval_trees

best_overall_eval_tree_file_paths      = snakemake.input.best_overall_eval_tree
iqtree_significance_summary_file_paths = snakemake.input.iqtree_significance_summary

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
        rfdist_log_file_path            = treesearch_rfDist_log_file_paths[i],
        best_tree_file_path             = best_treesearch_tree_file_paths[i],
        all_treesearch_trees_file_path  = treesearch_trees_file_paths[i],
        best_eval_tree_file_path        = best_eval_tree_log_file_paths[i],
        command                         = snakemake.params.raxml_command,
        all_eval_trees_file_path        = eval_trees_file_paths[i],
        rfdistances_file_path           = treesearch_rfDist_file_paths[i],
    )
    # fmt: on
    raxml_objects.append(raxml)

    # fmt:off
    raxml_db = db.Raxmlng.create(
        blmin   = raxml.blmin,
        blmax   = raxml.blmax,
        lh_eps  = raxml.lh_eps,
        raxml_param_epsilon     = raxml.model_param_epsilon,
        branch_length_smoothing = raxml.branch_length_smoothing,
        num_pars_trees          = raxml.num_pars_trees,
        num_rand_trees          = raxml.num_rand_trees,
        best_treesearch_llh     = raxml.best_treesearch_llh,
        best_evaluation_llh     = raxml.best_evaluation_llh,
        treesearch_total_time   = raxml.treesearch_total_time,
        avg_abs_rfdist_treesearch   = raxml.avg_abs_rfdist_treesearch,
        avg_rel_rfdist_treesearch   = raxml.avg_rel_rfdist_treesearch,
        num_unique_topos_treesearch = raxml.num_unique_topos_treesearch,
    )
    # fmt: on

    # RaxmlTreesearchTree
    raxml_treesearch_tree_db_objects = []
    raxml.db_best_treesearch_tree_object = None

    for tree_idx in range(raxml.get_num_of_trees()):
        tree_values = {}
        tree_values["llh"] = raxml.get_treesearch_llh_for_tree_index(tree_idx)
        tree_values["compute_time"] = raxml.get_treesearch_compute_time_for_tree_index(
            tree_idx
        )
        tree_values["newick_tree"] = raxml.get_newick_tree_for_tree_index(tree_idx)

        is_best = (
            raxml.tree_for_index_is_best(tree_idx)
            and not raxml.db_best_treesearch_tree_object
        )
        tree_values["is_best"] = is_best
        tree_values["number_of_taxa"] = raxml.get_number_of_taxa_for_tree_index(tree_idx)
        tree_values["total_branch_length"] = raxml.get_total_branch_length_for_tree_index(tree_idx)
        tree_values["average_branch_length"] = raxml.get_average_branch_length_for_tree_index(tree_idx)

        tree_values["program"] = raxml_db
        tree_values["seed"] = raxml.get_treesearch_seed_for_tree_index(tree_idx)

        raxml_treesearch_tree = db.RaxmlTreesearchTree.create(**tree_values)
        raxml_treesearch_tree_db_objects.append(raxml_treesearch_tree)

        if is_best:
            raxml.db_best_treesearch_tree_object = raxml_treesearch_tree

        # RFDistTreesearchTree
        # we need to create rfdistance object for each pairs of trees
        # at this point we can reference all trees from 0 to i
        # => we can create rfdistances for tree i with respect to all trees < i
        insert_into_rfdistance = []

        for tree_idx2 in range(tree_idx):
            rfdist_values = {}
            rfdist_values["tree1"] = raxml_treesearch_tree_db_objects[tree_idx2]
            rfdist_values["tree2"] = raxml_treesearch_tree_db_objects[tree_idx]
            rfdist_values["plain_rfdist"] = raxml.get_plain_rfdist_for_trees(
                (tree_idx2, tree_idx)
            )
            rfdist_values["normalized_rfdist"] = raxml.get_normalized_rfdist_for_trees(
                (tree_idx2, tree_idx)
            )

            insert_into_rfdistance.append(rfdist_values)

        with db.raxml_db.atomic():
            for batch in chunked(insert_into_rfdistance, 100):
                db.RFDistTreesearchTree.insert_many(batch).execute()

    # RaxmlEvalTree for best RaxmlTreesearchTree (raxml.db_best_treesearch_tree_object)
    for eval_tree_idx in range(raxml.get_num_of_eval_trees()):
        is_best = raxml.eval_tree_for_index_is_best(eval_tree_idx)
        # fmt: off
        eval_tree_values = {}
        eval_tree_values["start_tree"]      = raxml.db_best_treesearch_tree_object
        eval_tree_values["llh"]             = raxml.get_eval_llh_for_tree_index(eval_tree_idx)
        eval_tree_values["newick_tree"]     = raxml.get_newick_eval_tree_for_tree_index(eval_tree_idx)
        eval_tree_values["compute_time"]    = raxml.get_eval_compute_time_for_tree_index(eval_tree_idx)
        eval_tree_values["is_best"]         = is_best
        eval_tree_values["number_of_taxa"]          = raxml.get_number_of_taxa_for_eval_tree_index(eval_tree_idx)
        eval_tree_values["total_branch_length"]     = raxml.get_total_branch_length_for_eval_tree_index(eval_tree_idx)
        eval_tree_values["average_branch_length"]   = raxml.get_average_branch_length_for_eval_tree_index(eval_tree_idx)
        eval_tree_values["eval_blmin"]      = raxml.get_eval_blmin_for_tree_index(eval_tree_idx)
        eval_tree_values["eval_blmax"]      = raxml.get_eval_blmax_for_tree_index(eval_tree_idx)
        eval_tree_values["eval_lh_eps"]      = raxml.get_eval_lh_eps_for_tree_index(eval_tree_idx)
        eval_tree_values["eval_raxml_param_epsilon"]      = raxml.get_eval_model_param_epsilon_for_tree_index(eval_tree_idx)
        eval_tree_values["eval_raxml_brlen_smoothings"]      = raxml.get_eval_raxml_brlen_smoothings_for_tree_index(eval_tree_idx)
        # fmt: on
        raxml_eval_tree = db.RaxmlEvalTree.create(**eval_tree_values)

        if is_best:
            raxml.db_best_eval_tree = raxml_eval_tree
            best_eval_tree_objects.append(raxml_eval_tree)

# RFDistEvalTree
# fmt: off
experiment = create_Experiment(
    raxml_best_trees_path               = all_best_treesearch_trees_file_paths,
    raxml_best_eval_trees_path          = all_best_eval_trees_file_paths,
    rfdist_raxml_best_trees_path        = rfDist_best_treesearch_trees_file_paths,
    rfdist_raxml_best_eval_trees_path   = rfDist_best_eval_trees_file_paths,
    best_overall_eval_tree_file_path    = best_overall_eval_tree_file_paths,
    iqtree_statstest_results_file_path  = iqtree_significance_summary_file_paths,
)
# fmt: on

insert_into_rf_evaldist = []
for tree_idx1 in range(num_runs):
    raxml1 = raxml_objects[tree_idx1]
    tree1 = raxml1.db_best_eval_tree

    for tree_idx2 in range(tree_idx1 + 1, num_runs):
        raxml2 = raxml_objects[tree_idx2]
        # fmt: off
        rf_dist_values = {}
        rf_dist_values["tree1"] = tree1
        rf_dist_values["tree2"] = raxml2.db_best_eval_tree
        rf_dist_values["plain_rfdist"]      = experiment.get_plain_rfdist_for_raxml_eval_trees((tree_idx1, tree_idx2))
        rf_dist_values["normalized_rfdist"] = experiment.get_normalized_rfdist_for_raxml_eval_trees((tree_idx1, tree_idx2))
        # fmt: on

        insert_into_rf_evaldist.append(rf_dist_values)


with db.raxml_db.atomic():
    for batch in chunked(insert_into_rf_evaldist, 100):
        db.RFDistEvalTree.insert_many(batch).execute()

# Iqtree significance tests
best_overall_eval_tree = [tree for tree in best_eval_tree_objects if experiment.eval_tree_is_overall_best(tree.newick_tree)]
assert len(best_overall_eval_tree) >= 1, "Overall best eval tree not in all eval trees."
best_overall_eval_tree = best_overall_eval_tree[0]

insert_into_significance_table = []

for eval_tree in best_eval_tree_objects:
    newick_str = eval_tree.newick_tree
    statstest_values = {}
    statstest_values["reference_tree_id"] = best_overall_eval_tree
    statstest_values["tree_id"] = eval_tree
    statstest_values["iqtree_llh"] = experiment.get_iqtree_llh_for_eval_tree(newick_str)
    statstest_values["deltaL"] = experiment.get_iqtree_deltaL_for_eval_tree(newick_str)

    tests = experiment.get_iqtree_test_results_for_eval_tree(newick_str)

    if "bp-RELL" in tests:
        statstest_values["bpRell"] = tests["bp-RELL"]["score"]
        statstest_values["bpRell_significant"] = tests["bp-RELL"]["significant"]

    if "p-KH" in tests:
        statstest_values["pKH"] = tests["p-KH"]["score"]
        statstest_values["pKH_significant"] = tests["p-KH"]["significant"]

    if "p-SH" in tests:
        statstest_values["pSH"] = tests["p-SH"]["score"]
        statstest_values["pSH_significant"] = tests["p-SH"]["significant"]

    if "p-WKH" in tests:
        statstest_values["pWKH"] = tests["p-WKH"]["score"]
        statstest_values["pWKH_significant"] = tests["p-WKH"]["significant"]

    if "p-WSH" in tests:
        statstest_values["pWSH"] = tests["p-WSH"]["score"]
        statstest_values["pWSH_significant"] = tests["p-WSH"]["significant"]

    if "c-ELW" in tests:
        statstest_values["cELW"] = tests["c-ELW"]["score"]
        statstest_values["cELW_significant"] = tests["c-ELW"]["significant"]

    if "p-AU" in tests:
        statstest_values["pAU"] = tests["p-AU"]["score"]
        statstest_values["pAU_significant"] = tests["p-AU"]["significant"]

    insert_into_significance_table.append(statstest_values)

with db.raxml_db.atomic():
    for batch in chunked(insert_into_significance_table, 100):
        db.RaxmlEvalTreeStatsTest.insert_many(batch).execute()