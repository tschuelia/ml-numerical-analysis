import database as db
from experiment_parser import create_Run, create_Experiment
from raxml_parser import create_raxml
from iqtree_parser import create_iqtree

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

iqtree_treesearch_log_file_paths        = snakemake.input.iqtree_treesearch_log
iqtree_best_treesearch_tree_file_paths  = snakemake.input.iqtree_best_treesearch_tree
iqtree_treesearch_trees_file_paths      = snakemake.input.iqtree_treesearch_trees

iqtree_eval_log_file_paths          = snakemake.input.iqtree_eval_log
iqtree_best_eval_tree_file_paths    = snakemake.input.iqtree_best_eval_tree
iqtree_eval_trees_file_paths        = snakemake.input.iqtree_eval_trees

raxml_all_best_treesearch_trees_file_paths      = snakemake.input.raxml_all_best_treesearch_trees[0]
rfDist_raxml_best_treesearch_trees_file_paths   = snakemake.input.rfDist_raxml_best_treesearch_trees[0]

raxml_all_best_eval_trees_file_paths            = snakemake.input.raxml_all_best_eval_trees[0]
rfDist_raxml_best_eval_trees_file_paths         = snakemake.input.rfDist_raxml_best_eval_trees[0]

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
    run.raxmlng = raxmlng
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
    raxmlng.db_best_treesearch_tree_object = None
    for tree_idx in range(raxmlng.get_num_of_trees()):
        tree_values = {}
        tree_values["llh"] = raxmlng.get_treesearch_llh_for_tree_index(tree_idx)
        tree_values[
            "compute_time"
        ] = raxmlng.get_treesearch_compute_time_for_tree_index(tree_idx)
        tree_values["newick_tree"] = raxmlng.get_newick_tree_for_tree_index(tree_idx)

        is_best = (
            raxmlng.tree_for_index_is_best(tree_idx)
            and not raxmlng.db_best_treesearch_tree_object
        )
        tree_values["is_best"] = is_best

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

        if is_best:
            raxmlng.db_best_treesearch_tree_object = raxml_treesearch_tree

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

    # RaxmlEvalTree for best RaxmlTreesearchTree (raxmlng.db_best_treesearch_tree_object)
    for eval_tree_idx in range(raxmlng.get_num_of_eval_trees()):
        # fmt: off
        eval_tree_values = {}
        eval_tree_values["start_tree"]      = raxmlng.db_best_treesearch_tree_object
        eval_tree_values["llh"]             = raxmlng.get_eval_llh_for_tree_index(eval_tree_idx)
        eval_tree_values["newick_tree"]     = raxmlng.get_newick_eval_tree_for_tree_index(eval_tree_idx)
        eval_tree_values["compute_time"]    = raxmlng.get_eval_compute_time_for_tree_index( eval_tree_idx)
        eval_tree_values["is_best"]         = raxmlng.eval_tree_for_index_is_best(eval_tree_idx)
        eval_tree_values["eval_blmin"]      = raxmlng.get_eval_blmin_for_tree_index(eval_tree_idx)
        eval_tree_values["eval_blmax"]      = raxmlng.get_eval_blmax_for_tree_index(eval_tree_idx)
        # fmt: on
        raxml_eval_tree = db.RaxmlEvalTree.create(**eval_tree_values)

        if is_best:
            raxmlng.db_best_eval_tree = raxml_eval_tree


# RFDistEvalTree
# fmt: off
experiment = create_Experiment(
    raxml_best_trees_path               =raxml_all_best_treesearch_trees_file_paths,
    raxml_best_eval_trees_path          =raxml_all_best_eval_trees_file_paths,
    rfdist_raxml_best_trees_path        =rfDist_raxml_best_treesearch_trees_file_paths,
    rfdist_raxml_best_eval_trees_path   =rfDist_raxml_best_eval_trees_file_paths,
)
# fmt: on

insert_into_rf_evaldist = []
for tree_idx1 in range(len(run_python_objects)):
    raxmlng1 = run_python_objects[tree_idx1].raxmlng
    tree1 = raxmlng1.db_best_eval_tree

    for tree_idx2 in range(tree_idx1 + 1, len(run_python_objects)):
        raxmlng2 = run_python_objects[tree_idx2].raxmlng
        # fmt: off
        rf_dist_values = {}
        rf_dist_values["tree1"] = tree1
        rf_dist_values["tree2"] = raxmlng2.db_best_eval_tree
        rf_dist_values["plain_rfdist"]      = experiment.get_plain_rfdist_for_raxml_eval_trees((tree_idx1, tree_idx2))
        rf_dist_values["normalized_rfdist"] = experiment.get_normalized_rfdist_for_raxml_eval_trees((tree_idx1, tree_idx2))
        # fmt: on

        insert_into_rf_evaldist.append(rf_dist_values)


with db.db.atomic():
    db.RFDistEvalTree.insert_many(insert_into_rf_evaldist).execute()


#################################
# create Iqtree related
#################################
for run in run_python_objects:

    # Iqtree Program
    # fmt: off
    iqtree = create_iqtree(
        parameter_file_path             = params_file_paths[i],
        treesearch_log_file_path        = iqtree_treesearch_log_file_paths[i],
        eval_log_file_path              = iqtree_eval_log_file_paths[i],
        best_tree_file_path             = iqtree_best_treesearch_tree_file_paths[i],
        all_treesearch_trees_file_path  = iqtree_treesearch_trees_file_paths[i],
        best_eval_tree_file_path        = iqtree_best_eval_tree_file_paths[i],
        all_eval_trees_file_path        = iqtree_eval_trees_file_paths[i],
    )
   

    run.db_iqtree_object = db.Iqtree.create(
        run                     = run.db_run_object,
        num_pars_trees          = iqtree.num_pars_trees,
        num_rand_trees          = iqtree.num_rand_trees,
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

        tree_values["program"]  = run.db_iqtree_object
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
        eval_tree_values["eval_blmin"]  = iqtree.get_eval_blmin_for_tree_index(eval_tree_idx)
        eval_tree_values["eval_blmax"]  = iqtree.get_eval_blmax_for_tree_index(eval_tree_idx)
        # fmt: on
        iqtree_eval_tree = db.IqtreeEvalTree.create(**eval_tree_values)
