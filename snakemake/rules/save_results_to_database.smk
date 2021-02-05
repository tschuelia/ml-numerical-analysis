rule save_results_to_database:
    input:
        params_file             = expand(f"{full_dir}parameters.json", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
        
        #raxml_treesearch_log        = expand(f"{full_file_path_raxml}.allTreesearchLogs", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
        #raxml_best_treesearch_tree  = expand(f"{full_file_path_raxml}.bestTreeOfRun", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
        #raxml_treesearch_trees      = expand(f"{full_file_path_raxml}.allTreesCollected", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
        
        raxml_eval_log          = expand(f"{full_file_path_raxml}.allEvalLogs", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
        raxml_best_eval_tree    = expand(f"{full_file_path_raxml}.bestEvalTreeOfRun", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
        raxml_eval_trees        = expand(f"{full_file_path_raxml}.allEvalTreesCollected", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),

        raxml_treesearch_rfDist_log = expand(f"{full_file_path_raxml}.raxml.rfDistances.log", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
        raxml_treesearch_rfDist     = expand(f"{full_file_path_raxml}.raxml.rfDistances", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),

        raxml_iqtree_statstest_results = expand(f"{full_file_path_iqtree}.iqtree", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),

        iqtree_treesearch_log       = expand(f"{full_file_path_iqtree}.allTreesearchLogs", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
        iqtree_best_treesearch_tree = expand(f"{full_file_path_iqtree}.bestTreeOfRun", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
        iqtree_treesearch_trees     = expand(f"{full_file_path_iqtree}.allTreesCollected", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),

        iqtree_eval_log         = expand(f"{full_file_path_iqtree}.allEvalLogs", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
        iqtree_best_eval_tree   = expand(f"{full_file_path_iqtree}.bestEvalTreeOfRun", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
        iqtree_eval_trees       = expand(f"{full_file_path_iqtree}.allEvalTreesCollected", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),

        raxml_all_best_treesearch_trees     = expand(f"{outdir}/bestTreesCollected", outdir=outdir),
        rfDist_raxml_best_treesearch_trees  = expand(f"{outdir}/bestTrees.raxml.rfDistances", outdir=outdir),
        raxml_all_best_eval_trees           = expand(f"{outdir}/bestEvalTreesCollected", outdir=outdir),
        rfDist_raxml_best_eval_trees        = expand(f"{outdir}/bestEvalTrees.raxml.rfDistances", outdir=outdir),

    output:
        database = f"{outdir}/results.sqlite3", 
    params:
        raxml_command = config["parameters"]["raxml-ng"]["command"],
    script:
        "scripts/save_results_to_database.py"

