rule save_fasttree_results_to_database:
    input:
        params_file  = expand(f"{full_dir_fasttree}parameters.json", blmin=blmin_opts, lh_eps=lh_eps_opts),

        ##### treesearch
        treesearch_trees    = expand(f"{full_file_path_fasttree}.allTreesCollected", blmin=blmin_opts, lh_eps=lh_eps_opts),
        treesearch_log      = expand(f"{full_file_path_fasttree}.allTreesearchLogs", blmin=blmin_opts, lh_eps=lh_eps_opts),
        best_treesearch_tree = expand(f"{full_file_path_fasttree}.bestTreeOfRun", blmin=blmin_opts, lh_eps=lh_eps_opts),

        ##### iqtree significance tests
        filtered_trees_clusters     = f"{base_dir_fasttree}filteredTreesClusters",
        iqtree_statstests_results   =f"{base_dir_fasttree}significance.iqtree",
    output:
        database = f"{base_dir}fasttree_results.sqlite3"
    script:
        "scripts/save_fasttree_results_to_database.py"
