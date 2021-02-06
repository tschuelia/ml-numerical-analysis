rule save_raxml_results_to_database:
    input:
        params_file      = expand(f"{full_dir_raxml}parameters.json", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, raxml_param_epsilon=raxml_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts),
        ##### treesearch
        # all infered trees for one combination of parameters collected in a single file
        treesearch_trees = expand(f"{full_file_path_raxml}.allTreesCollected", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, raxml_param_epsilon=raxml_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts),

        # all logs for one run
        treesearch_log = expand(f"{full_file_path_raxml}.allTreesearchLogs", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, raxml_param_epsilon=raxml_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts),

        # best tree per run
        best_treesearch_tree = expand(f"{full_file_path_raxml}.bestTreeOfRun", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, raxml_param_epsilon=raxml_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts),
       
        ##### evaluation
        _               = expand(f"{full_file_path_raxml_eval}.raxml.eval.log", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, raxml_param_epsilon=raxml_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts, 
                                                              blmin_eval=blmin_opts, blmax_eval=blmax_opts, lh_eps_eval=lh_eps_opts, raxml_param_epsilon_eval=raxml_param_epsilon_opts, raxml_brlen_smoothings_eval=raxml_brlen_smoothings_opts),
        eval_trees      = expand(f"{full_file_path_raxml}.allEvalTreesCollected", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, raxml_param_epsilon=raxml_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts),
        eval_log        = expand(f"{full_file_path_raxml}.allEvalLogs", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, raxml_param_epsilon=raxml_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts),
        best_eval_tree  = expand(f"{full_file_path_raxml}.bestEvalTreeOfRun", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, raxml_param_epsilon=raxml_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts),

        ##### rfdistances
        treesearch_rfDist_log   = expand(f"{full_file_path_raxml}.raxml.rfDistances.log", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, raxml_param_epsilon=raxml_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts),
        treesearch_rfDist       = expand(f"{full_file_path_raxml}.raxml.rfDistances", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, raxml_param_epsilon=raxml_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts),

        all_best_treesearch_trees     = f"{base_dir_raxml}bestTreesCollected",
        rfDist_best_treesearch_trees  = f"{base_dir_raxml}bestTrees.raxml.rfDistances",

        all_best_eval_trees           = f"{base_dir_raxml}bestEvalTreesCollected",
        rfDist_best_eval_trees        = f"{base_dir_raxml}bestEvalTrees.raxml.rfDistances",
    output:
        database = f"{base_dir}raxml_results.sqlite3", 
    params:
        raxml_command = config["parameters"]["raxml-ng"]["command"],
    script:
        "scripts/save_raxml_results_to_database.py"