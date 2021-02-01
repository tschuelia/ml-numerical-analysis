import database as db
from experiment_parser import create_Run, create_Experiment
from raxml_parser import create_raxml

# initialize empty database
db.db.init(snakemake.output.database)
db.db.connect()
db.db.create_tables(
    [
        db.Run,
        db.Raxmlng,
        db.Iqtree,
        db.RaxmlTreesearchTree,
        db.IqtreeTreesearchTree,
        db.RaxmlEvalTree,
        db.IqtreeEvalTree,
        db.RFDistTreesearchTree,
        db.RFDistEvalTree,
    ]
)

# fmt: off
params_file_paths = snakemake.input.params_file

raxml_treesearch_log_file_paths         = snakemake.input.raxml_treesearch_log
raxml_best_treesearch_tree_file_paths   = snakemake.input.raxml_best_treesearch_tree
raxml_treesearch_trees_file_paths       = snakemake.input.raxml_treesearch_trees

raxml_eval_log_file_paths           = snakemake.input.raxml_eval_log
raxml_best_eval_tree_log_file_paths = snakemake.input.raxml_best_eval_tree
raxml_eval_trees_file_paths         = snakemake.input.raxml_eval_trees

raxml_iqtree_statstest_results_file_paths = snakemake.input.raxml_iqtree_statstest_results

raxml_treesearch_rfDist_log_file_paths  = snakemake.input.raxml_treesearch_rfDist_log
raxml_treesearch_rfDist_file_paths      = snakemake.input.raxml_treesearch_rfDist
# fmt: on

num_runs = len(params_file_paths)
#################################
# create Run objects
#################################
run_python_objects = []

for i in range(num_runs):
    run_python_objects.append(create_Run(params_file_paths[i]))

for run in run_python_objects:
    run.db_run_object = db.Run.create(blmin=run.blmin, blmax=run.blmax)

#################################
# create Raxmlng related
#################################

for run in run_python_objects:
    # fmt:off
    raxmlng = create_raxml(
        parameter_file_path             = params_file_paths[i],
        treesearch_log_file_path        = raxml_treesearch_log_file_paths[i],
        eval_log_file_path              = raxml_eval_log_file_paths[i],
        rfdist_log_file_path            = raxml_treesearch_rfDist_log_file_paths[i],
        best_tree_file_path             = raxml_best_treesearch_tree_file_paths[i],
        all_treesearch_trees_file_path  = raxml_treesearch_trees_file_paths[i],
        iqtree_statstest_results_file_path  = raxml_iqtree_statstest_results_file_paths[i],
        best_eval_tree_file_path        = raxml_best_eval_tree_log_file_paths[i],
        command                         = snakemake.params.raxml_command,
        all_eval_trees_file_path        = raxml_eval_trees_file_paths[i],
        rfdistances_file_path           = raxml_treesearch_rfDist_file_paths[i],
    )
    # fmt: on

    # Raxmlng Program
    run.db_raxml_object = db.Raxmlng.create(
        run=run.db_run_object,
        num_pars_trees=raxmlng.num_pars_trees,
        num_rand_trees=raxmlng.num_rand_trees,
        best_treesearch_llh=raxmlng.best_treesearch_llh,
        best_evaluation_llh=raxmlng.best_evaluation_llh,
        treesearch_total_time=raxmlng.treesearch_total_time,
        avg_abs_rfdist_treesearch=raxmlng.avg_abs_rfdist_treesearch,
        avg_rel_rfdist_treesearch=raxmlng.avg_rel_rfdist_treesearch,
        num_unique_topos_treesearch=raxmlng.num_unique_topos_treesearch,
    )

    # RaxmlTreesearchTree
    raxml_treesearch_tree_db_objects = []

    for tree_idx in range(raxmlng.get_num_of_tres()):
        tree_values = {}
        tree_values["llh"] = raxmlng.get_treesearch_llh_for_tree_index(tree_idx)
        tree_values[
            "compute_time"
        ] = raxmlng.get_treesearch_compute_time_for_tree_index(tree_idx)
        tree_values["newick_tree"] = raxmlng.get_newick_tree_for_tree_index(tree_idx)
        tree_values["is_best"] = raxmlng.tree_for_index_is_best(tree_idx)

        tree_values["program"] = run.db_raxml_object
        tree_values["seed"] = raxmlng.get_treesearch_seed_for_tree_index(tree_idx)

        tree_values["iqtree_llh"] = raxmlng.get_iqtree_llh_for_tree_index(tree_idx)

        tree_values["deltaL"] = raxmlng.get_iqtree_deltaL_for_tree_index(tree_idx)

        tests = raxmlng.get_iqtree_test_results_for_tree_index(tree_idx)

        if "bp-RELL" in tests:
            tree_values["bpRell"] = tests["bp-RELL"]["score"]
            tree_values["bpRell_significant"] = tests["bp-RELL"]["significant"]

        if "p-KH" in tests:
            tree_values["pKH"] = tests["p-KH"]["score"]
            tree_values["pKH_significant"] = tests["p-KH"]["significant"]

        if "p-SH" in tests:
            tree_values["pSH"] = tests["p-SH"]["score"]
            tree_values["pSH_significant"] = tests["p-SH"]["significant"]

        if "p-WKH" in tests:
            tree_values["pWKH"] = tests["p-WKH"]["score"]
            tree_values["pWKH_significant"] = tests["p-WKH"]["significant"]

        if "p-WSH" in tests:
            tree_values["pWSH"] = tests["p-WSH"]["score"]
            tree_values["pWSH_significant"] = tests["p-WSH"]["significant"]

        if "c-ELW" in tests:
            tree_values["cELW"] = tests["c-ELW"]["score"]
            tree_values["cELW_significant"] = tests["c-ELW"]["significant"]

        if "p-AU" in tests:
            tree_values["pAU"] = tests["p-AU"]["score"]
            tree_values["pAU_significant"] = tests["p-AU"]["significant"]

        raxml_treesearch_tree = db.RaxmlTreesearchTree.create(**tree_values)
        raxml_treesearch_tree_db_objects.append(raxml_treesearch_tree)

        # RFDistTreesearchTree
        # we need to create rfdistance object for each pairs of trees
        # at this point we can reference all trees from 0 to i
        # => we can create rfdistances for tree i with respect to all trees < i
        insert_into_rfdistance = []

        for tree_idx2 in range(tree_idx):
            rfdist_values = {}
            rfdist_values["tree1"] = raxml_treesearch_tree_db_objects[tree_idx2]
            rfdist_values["tree2"] = raxml_treesearch_tree_db_objects[tree_idx]
            rfdist_values["plain_rfdist"] = raxmlng.get_plain_rfdist_for_trees(
                (tree_idx2, tree_idx)
            )
            rfdist_values[
                "normalized_rfdist"
            ] = raxmlng.get_normalized_rfdist_for_trees((tree_idx2, tree_idx))

            insert_into_rfdistance.append(rfdist_values)

        with db.db.atomic():
            db.RFDistTreesearchTree.insert_many(insert_into_rfdistance).execute()
