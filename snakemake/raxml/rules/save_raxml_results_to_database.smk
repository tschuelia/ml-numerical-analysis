rule save_raxml_results_to_database:
    input:
        params_file      = expand(f"{full_dir_raxml}parameters.json", lheps=lheps_opts, lheps_trip=lheps_trip_opts),
        ##### treesearch
        # all infered trees for one combination of parameters collected in a single file
        treesearch_trees    = expand(f"{full_file_path_raxml}.allTreesCollected", lheps=lheps_opts, lheps_trip=lheps_trip_opts),

        # all logs for one run
        treesearch_log = expand(f"{full_file_path_raxml}.allTreesearchLogs", lheps=lheps_opts, lheps_trip=lheps_trip_opts),

        # all starting tree files for one run
        starting_trees      = expand(f"{full_file_path_raxml}.allStartingTreesCollected", lheps=lheps_opts, lheps_trip=lheps_trip_opts),
        starting_eval_trees = expand(f"{full_file_path_raxml}.allStartingEvalTreesCollected", lheps=lheps_opts, lheps_trip=lheps_trip_opts),
        starting_eval_logs  = expand(f"{full_file_path_raxml}.allStartingEvalLogsCollected", lheps=lheps_opts, lheps_trip=lheps_trip_opts),

        # best tree per run
        best_treesearch_tree = expand(f"{full_file_path_raxml}.bestTreeOfRun", lheps=lheps_opts, lheps_trip=lheps_trip_opts),

        ##### evaluation
        _               = expand(f"{full_file_path_raxml_eval_pars}.raxml.eval.log", seed=pars_seeds, lheps=lheps_opts, lheps_trip=lheps_trip_opts),
        __              = expand(f"{full_file_path_raxml_eval_rand}.raxml.eval.log", seed=rand_seeds, lheps=lheps_opts, lheps_trip=lheps_trip_opts),

        eval_trees      = expand(f"{full_file_path_raxml}.allEvalTreesCollected", lheps=lheps_opts, lheps_trip=lheps_trip_opts),
        eval_log        = expand(f"{full_file_path_raxml}.allEvalLogs", lheps=lheps_opts, lheps_trip=lheps_trip_opts),
        best_eval_tree  = expand(f"{full_file_path_raxml}.bestEvalTreeOfRun", lheps=lheps_opts, lheps_trip=lheps_trip_opts),

        ##### iqtree significance tests
        filtered_trees_clusters     = f"{base_dir_raxml}filteredEvalTreesClusters",
        iqtree_statstests_results   = f"{base_dir_raxml}significance.iqtree",

        # pairwise results
        pairwise_iqtree_statstests_results= f"{base_dir_raxml}pairwiseSignificanceEval.json",

        ##### iqtree significance tests for eval and starting
        filtered_trees_clusters_ev_starting     = f"{base_dir_raxml}filteredEvalAndStartingTreesClusters",
        iqtree_statstests_results_ev_starting   = f"{base_dir_raxml}significanceEvalAndStarting.iqtree",

        # pairwise results
        pairwise_iqtree_statstests_results_ev_starting= f"{base_dir_raxml}pairwiseSignificanceEvalAndStarting.json",
    output:
        database = f"{base_dir}raxml_results.sqlite3", 
    params:
        raxml_command = config["parameters"]["raxml-ng"]["command"],
    script:
        "scripts/save_raxml_results_to_database.py"