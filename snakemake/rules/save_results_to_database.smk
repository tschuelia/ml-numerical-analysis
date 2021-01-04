rule save_results_to_database:
    input:
        params_file             = expand(f"{full_dir}parameters.json", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
        best_tree_raxml         = expand(f"{full_file_path_raxml}.raxml.bestTree", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
        all_trees_raxml         = expand(f"{full_file_path_raxml}.raxml.mlTrees", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
        iqtree_results          = expand(f"{full_file_path_iqtree}.iqtree", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
        iqtree_trees            = expand(f"{full_file_path_iqtree}.trees", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
        raxml_treesearch_log    = expand(f"{full_file_path_raxml}.raxml.treesearch.log", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
        iqtree_test_log         = expand(f"{full_file_path_iqtree}.iqtree_tests.log", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
        rfDistances_log         = expand(f"{full_file_path_raxml}.raxml.rfDistances.log", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
        rfDistances             = expand(f"{full_file_path_raxml}.raxml.rfDistances", blmin=blmin_opts, blmax=blmax_opts, outdir=outdir),
    output:
        database = f"{outdir}/results.sqlite3",    
    script:
        "scripts/save_results_to_database.py"