rule save_raxml_results_to_database:
    input:
        params_file      = expand(f"{full_dir_raxml}parameters.json", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, model_epsilon=model_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts, spr_lh_eps=spr_lh_eps_opts, bfgs_fac=bfgs_fac_opts,),
        ##### treesearch
        # all infered trees for one combination of parameters collected in a single file
        treesearch_trees = expand(f"{full_file_path_raxml}.allTreesCollected", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, model_epsilon=model_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts, spr_lh_eps=spr_lh_eps_opts, bfgs_fac=bfgs_fac_opts,),

        # all logs for one run
        treesearch_log = expand(f"{full_file_path_raxml}.allTreesearchLogs", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, model_epsilon=model_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts, spr_lh_eps=spr_lh_eps_opts, bfgs_fac=bfgs_fac_opts,),

        # best tree per run
        best_treesearch_tree = expand(f"{full_file_path_raxml}.bestTreeOfRun", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, model_epsilon=model_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts, spr_lh_eps=spr_lh_eps_opts, bfgs_fac=bfgs_fac_opts,),
       
        ##### evaluation
        _               = expand(f"{full_file_path_raxml_eval_pars}.raxml.eval.log", seed=pars_seeds, blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, model_epsilon=model_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts, spr_lh_eps=spr_lh_eps_opts, bfgs_fac=bfgs_fac_opts),
        __              = expand(f"{full_file_path_raxml_eval_rand}.raxml.eval.log", seed=rand_seeds, blmin=blmin_opts,blmax=blmax_opts,lh_eps=lh_eps_opts,model_epsilon=model_param_epsilon_opts,raxml_brlen_smoothings=raxml_brlen_smoothings_opts,spr_lh_eps=spr_lh_eps_opts,bfgs_fac=bfgs_fac_opts),

        eval_trees      = expand(f"{full_file_path_raxml}.allEvalTreesCollected", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, model_epsilon=model_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts, spr_lh_eps=spr_lh_eps_opts, bfgs_fac=bfgs_fac_opts,),
        eval_log        = expand(f"{full_file_path_raxml}.allEvalLogs", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, model_epsilon=model_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts, spr_lh_eps=spr_lh_eps_opts, bfgs_fac=bfgs_fac_opts,),
        best_eval_tree  = expand(f"{full_file_path_raxml}.bestEvalTreeOfRun", blmin=blmin_opts, blmax=blmax_opts, lh_eps=lh_eps_opts, model_epsilon=model_param_epsilon_opts, raxml_brlen_smoothings=raxml_brlen_smoothings_opts, spr_lh_eps=spr_lh_eps_opts, bfgs_fac=bfgs_fac_opts,),

        ##### iqtree significance tests
        filtered_trees_clusters     = f"{base_dir_raxml}filteredEvalTreesClusters",
        iqtree_statstests_results   = f"{base_dir_raxml}significance.iqtree",
    output:
        database = f"{base_dir}raxml_results.sqlite3", 
    params:
        raxml_command = config["parameters"]["raxml-ng"]["command"],
        blmin_eval= blmin_eval,
        blmax_eval=blmax_eval,
        lh_eps_eval=lh_eps_eval,
        model_param_epsilon_eval=model_param_epsilon_eval,
        raxml_brlen_smoothings_eval=raxml_brlen_smoothings_eval,
        bfgs_fac_eval=bfgs_fac_eval,
    script:
        "scripts/save_raxml_results_to_database.py"