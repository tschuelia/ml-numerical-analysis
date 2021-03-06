rule save_iqtree_results_to_database:
    input: 
        params_file  = expand(f"{full_dir_iqtree}parameters.json", blmin=blmin_opts, blmax=blmax_opts, model_epsilon=model_param_epsilon_opts, lh_eps=lh_eps_opts),

        ##### treesearch
        treesearch_trees = expand(f"{full_file_path_iqtree}.allTreesCollected", blmin=blmin_opts, blmax=blmax_opts, model_epsilon=model_param_epsilon_opts, lh_eps=lh_eps_opts),
        treesearch_log = expand(f"{full_file_path_iqtree}.allTreesearchLogs", blmin=blmin_opts, blmax=blmax_opts, model_epsilon=model_param_epsilon_opts, lh_eps=lh_eps_opts),
        best_treesearch_tree = expand(f"{full_file_path_iqtree}.bestTreeOfRun", blmin=blmin_opts, blmax=blmax_opts, model_epsilon=model_param_epsilon_opts, lh_eps=lh_eps_opts),

        ##### evaluation
        _               = expand(f"{full_file_path_iqtree_eval}.iqtree.eval.log", seed=pars_seeds, blmin=blmin_opts, blmax=blmax_opts, model_epsilon=model_param_epsilon_opts, lh_eps=lh_eps_opts,),
        eval_trees      = expand(f"{full_file_path_iqtree}.allEvalTreesCollected", blmin=blmin_opts, blmax=blmax_opts, model_epsilon=model_param_epsilon_opts, lh_eps=lh_eps_opts),
        eval_log        = expand(f"{full_file_path_iqtree}.allEvalLogs", blmin=blmin_opts, blmax=blmax_opts, model_epsilon=model_param_epsilon_opts, lh_eps=lh_eps_opts),
        best_eval_tree  = expand(f"{full_file_path_iqtree}.bestEvalTreeOfRun", blmin=blmin_opts, blmax=blmax_opts, model_epsilon=model_param_epsilon_opts, lh_eps=lh_eps_opts),

        ##### iqtree significance tests
        filtered_trees_clusters     = f"{base_dir_iqtree}filteredEvalTreesClusters",
        iqtree_statstests_results   = f"{base_dir_iqtree}significance.iqtree",
    params:
        blmin_eval = blmin_eval,
        blmax_eval = blmax_eval,
        lh_eps_eval = lh_eps_eval,
        model_param_epsilon_eval = model_param_epsilon_eval,
    output:
        database = f"{base_dir}iqtree_results.sqlite3"
    script:
        "scripts/save_iqtree_results_to_database.py"

