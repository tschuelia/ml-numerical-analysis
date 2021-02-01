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
# Raxmlng Program

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