rule save_iqtree_results_to_database:
    input: 
        params_file  = expand(f"{full_dir_iqtree}parameters.json", lh_eps=lh_eps_opts),

        ##### treesearch
        treesearch_trees = expand(f"{full_file_path_iqtree}.allTreesCollected", lh_eps=lh_eps_opts),
        treesearch_log = expand(f"{full_file_path_iqtree}.allTreesearchLogs", lh_eps=lh_eps_opts),
        best_treesearch_tree = expand(f"{full_file_path_iqtree}.bestTreeOfRun", lh_eps=lh_eps_opts),

        ##### evaluation
        _               = expand(f"{full_file_path_iqtree_eval}.iqtree.eval.log", seed=pars_seeds, lh_eps=lh_eps_opts),
        eval_trees      = expand(f"{full_file_path_iqtree}.allEvalTreesCollected", lh_eps=lh_eps_opts),
        eval_log        = expand(f"{full_file_path_iqtree}.allEvalLogs", lh_eps=lh_eps_opts),
        best_eval_tree  = expand(f"{full_file_path_iqtree}.bestEvalTreeOfRun", lh_eps=lh_eps_opts),

        ##### iqtree significance tests
        filtered_trees_clusters     = f"{base_dir_iqtree}filteredEvalTreesClusters",
        iqtree_statstests_results   = f"{base_dir_iqtree}significance.iqtree",

        # pairwise results
        pairwise_iqtree_statstests_results=f"{base_dir_iqtree}pairwiseSignificanceEval.json",

    params:
    output:
        database = f"{base_dir}iqtree_results.sqlite3"
    script:
        "scripts/save_iqtree_results_to_database.py"

