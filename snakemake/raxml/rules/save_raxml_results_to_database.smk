rule save_raxml_results_to_database:
    input:
        params_file      = expand(f"{full_dir_raxml}parameters.json", lheps_auto=auto_opts, lheps_fast=fast_opts, lheps_slow=slow_opts, lheps_full=full_opts, lheps_trip=trip_opts),
        ##### treesearch
        # all infered trees for one combination of parameters collected in a single file
        treesearch_trees    = expand(f"{full_file_path_raxml}.allTreesCollected", lheps_auto=auto_opts, lheps_fast=fast_opts, lheps_slow=slow_opts, lheps_full=full_opts, lheps_trip=trip_opts),

        # all logs for one run
        treesearch_log = expand(f"{full_file_path_raxml}.allTreesearchLogs", lheps_auto=auto_opts, lheps_fast=fast_opts, lheps_slow=slow_opts, lheps_full=full_opts, lheps_trip=trip_opts),

        # all starting tree files for one run
        starting_trees      = expand(f"{full_file_path_raxml}.allStartingTreesCollected", lheps_auto=auto_opts, lheps_fast=fast_opts, lheps_slow=slow_opts, lheps_full=full_opts, lheps_trip=trip_opts),
        starting_eval_trees = expand(f"{full_file_path_raxml}.allStartingEvalTreesCollected", lheps_auto=auto_opts, lheps_fast=fast_opts, lheps_slow=slow_opts, lheps_full=full_opts, lheps_trip=trip_opts),
        starting_eval_logs  = expand(f"{full_file_path_raxml}.allStartingEvalLogsCollected", lheps_auto=auto_opts, lheps_fast=fast_opts, lheps_slow=slow_opts, lheps_full=full_opts, lheps_trip=trip_opts),

        # best tree per run
        best_treesearch_tree = expand(f"{full_file_path_raxml}.bestTreeOfRun", lheps_auto=auto_opts, lheps_fast=fast_opts, lheps_slow=slow_opts, lheps_full=full_opts, lheps_trip=trip_opts),

        ##### evaluation
        _               = expand(f"{full_file_path_raxml_eval_pars}.raxml.eval.log", seed=pars_seeds, lheps_auto=auto_opts, lheps_fast=fast_opts, lheps_slow=slow_opts, lheps_full=full_opts, lheps_trip=trip_opts),
        __              = expand(f"{full_file_path_raxml_eval_rand}.raxml.eval.log", seed=rand_seeds, lheps_auto=auto_opts, lheps_fast=fast_opts, lheps_slow=slow_opts, lheps_full=full_opts, lheps_trip=trip_opts),

        eval_trees      = expand(f"{full_file_path_raxml}.allEvalTreesCollected", lheps_auto=auto_opts, lheps_fast=fast_opts, lheps_slow=slow_opts, lheps_full=full_opts, lheps_trip=trip_opts),
        eval_log        = expand(f"{full_file_path_raxml}.allEvalLogs", lheps_auto=auto_opts, lheps_fast=fast_opts, lheps_slow=slow_opts, lheps_full=full_opts, lheps_trip=trip_opts),
        best_eval_tree  = expand(f"{full_file_path_raxml}.bestEvalTreeOfRun", lheps_auto=auto_opts, lheps_fast=fast_opts, lheps_slow=slow_opts, lheps_full=full_opts, lheps_trip=trip_opts),

        ##### iqtree significance tests
        filtered_trees_clusters     = f"{base_dir_raxml}filteredEvalTreesClusters",
        iqtree_statstests_results   = f"{base_dir_raxml}significance.iqtree",

        ##### iqtree significance tests for eval and starting
        filtered_trees_clusters_ev_starting     = f"{base_dir_raxml}filteredEvalAndStartingTreesClusters",
        iqtree_statstests_results_ev_starting   = f"{base_dir_raxml}significanceEvalAndStarting.iqtree",
    output:
        database = f"{base_dir}raxml_results.sqlite3", 
    params:
        raxml_command = config["parameters"]["raxml-ng"]["command"],
    script:
        "scripts/save_raxml_results_to_database.py"