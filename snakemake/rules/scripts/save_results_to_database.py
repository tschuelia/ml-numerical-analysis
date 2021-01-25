import database as db
from experiment_parser import create_Run, create_Experiment

# initialize empty database
db.db.init(snakemake.output.database)
db.db.connect()
db.db.create_tables([db.Run, db.Tree, db.RFDistance])

# snakemake.input.[something] is a list of filepaths
params_file_paths = snakemake.input.params_file
best_tree_file_paths = snakemake.input.best_tree_raxml
all_trees_raxml_file_paths = snakemake.input.all_trees_raxml
iqtree_results_file_paths = snakemake.input.iqtree_results
iqtree_trees_file_paths = snakemake.input.iqtree_trees
raxml_treesearch_log_file_paths = snakemake.input.raxml_treesearch_log
iqtree_test_log_file_paths = snakemake.input.iqtree_test_log
rfDistances_log_file_path = snakemake.input.rfDistances_log
rfDistance_path = snakemake.input.rfDistances

# these are single files only
rfDistances_best_trees = snakemake.input.rfDistances_best_trees[0]
best_trees_collected = snakemake.input.best_trees_collected[0]

num_runs = len(best_tree_file_paths)

run_python_objects = []

for i in range(num_runs):
    run_python_objects.append(
        create_Run(
            parameter_file_path=params_file_paths[i],
            best_raxml_tree_file_path=best_tree_file_paths[i],
            all_raxml_trees_file_path=all_trees_raxml_file_paths[i],
            raxml_treesearch_log_file_path=raxml_treesearch_log_file_paths[i],
            all_iqtree_trees_file_path=iqtree_trees_file_paths[i],
            iqtree_results_file_path=iqtree_results_file_paths[i],
            iqtree_test_log_file_path=iqtree_test_log_file_paths[i],
            rfdistances_file_path=rfDistance_path[i],
            raxml_rfdistance_logfile_path=rfDistances_log_file_path[i],
        )
    )

for run in run_python_objects:
    run.db_run_object = db.Run.create(
        num_raxml_pars_trees=run.get_num_raxml_pars_trees(),
        num_raxml_rand_trees=run.get_num_raxml_rand_trees(),
        blmin=run.get_blmin(),
        blmax=run.get_blmax(),
        average_absolute_rf_distance=run.get_average_absolute_rf_distance(),
        average_relative_rf_distance=run.get_average_relative_rf_distance(),
        num_unique_topos=run.get_num_unique_topos(),
        raxml_best_llh=run.get_best_raxml_llh(),
        iqtree_best_llh=run.get_best_iqtree_llh(),
        raxml_treesearch_elapsed_time=run.get_raxml_treesearch_elapsed_time(),
    )

    tree_objects = []

    for i in range(run.get_num_of_trees()):
        tree_values = {}
        tree_values["run"] = run.db_run_object
        tree_values["raxml_tree"] = run.get_raxml_tree_for_tree_index(i)
        tree_values["iqtree_tree"] = run.get_iqtree_tree_for_tree_index(i)
        tree_values["raxml_llh"] = run.get_raxml_llh_for_tree_index(i)
        tree_values["iqtree_llh"] = run.get_iqtree_llh_for_tree_index(i)
        tree_values["raxml_seed"] = run.get_raxml_seed_for_tree_index(i)

        is_best = run.tree_for_index_is_best(i)
        tree_values["is_best"] = is_best

        tree_values["deltaL"] = run.get_iqtree_deltaL_for_tree_index(i)

        tests = run.get_iqtree_test_results_for_tree_index(i)

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

        tree = db.Tree.create(**tree_values)
        tree_objects.append(tree)

        if is_best:
            run.db_best_tree = tree

        insert_into_rfdistance = []

        # we need to create rfdistance object for each pairs of trees
        # at this point we can reference all trees from 0 to i
        # => we can create rfdistances for tree i with respect to all trees < i
        for j in range(i):
            rf_dist_values = {}
            rf_dist_values["tree1"] = tree_objects[j]
            rf_dist_values["tree2"] = tree_objects[i]
            rf_dist_values["plain_rf_distance"] = run.get_plain_rfdistance_for_trees(
                (i, j)
            )
            rf_dist_values[
                "normalized_rf_distance"
            ] = run.get_normalized_rfdistance_for_trees((i, j))

            insert_into_rfdistance.append(rf_dist_values)

        with db.db.atomic():
            db.RFDistance.insert_many(insert_into_rfdistance).execute()

experiment = create_Experiment(
    runs=run_python_objects,
    best_trees_path=best_trees_collected,
    rfdist_best_trees_path=rfDistances_best_trees,
)

# create rfdistance db objects for best trees
insert_into_rfdistance = []
for i in range(len(run_python_objects)):
    tree1 = run_python_objects[i].db_best_tree
    for j in range(i + 1, len(run_python_objects)):
        rf_dist_values = {}
        rf_dist_values["tree1"] = tree1
        rf_dist_values["tree2"] = run_python_objects[j].db_best_tree
        rf_dist_values["plain_rf_distance"] = experiment.get_plain_rfdistance_for_trees(
            (i, j)
        )
        rf_dist_values[
            "normalized_rf_distance"
        ] = experiment.get_normalized_rfdistance_for_trees((i, j))

        insert_into_rfdistance.append(rf_dist_values)

with db.db.atomic():
    db.RFDistance.insert_many(insert_into_rfdistance).execute()
