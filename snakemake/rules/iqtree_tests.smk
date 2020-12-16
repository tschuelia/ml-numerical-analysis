rule perform_iqtree_tests:
    input: 
        msa         = config["data"]["input"],
        ml_trees    = f"{full_file_path_raxml}.raxml.mlTrees",
        best_tree   = f"{full_file_path_raxml}.raxml.bestTree",
    output:
        summary     = f"{full_file_path_iqtree}.iqtree",
        best_iqtree = f"{full_file_path_iqtree}.treefile",
        all_trees   = f"{full_file_path_iqtree}.trees",
    params:
        model       = config["parameters"]["model"]["iqtree"],
        threads     = config["parameters"]["iqtree"]["threads"],
        prefix      = full_file_path_iqtree,
    log:
        f"{full_file_path_iqtree}.iqtree_tests.log"
    shell:
        "{iqtree_command} "
        "-blmin {wildcards.blmin} "
        "-s {input.msa} "
        "-m {params.model} "
        "-pre {params.prefix} "
        "-z {input.ml_trees} "
        "-te {input.best_tree} "
        "-n 0 "
        "-zb 10000 "
        "-zw "
        "-au "
        "-nt {params.threads} "
        "> {log} "

rule filter_trees:
    """
    Construct tree set containing only trees that are not significantly
    worse than the best scoring tree according to IQ-Tree results.
    """
    input: 
        summary     = f"{full_file_path_iqtree}.iqtree",
        ml_trees    = f"{full_file_path_raxml}.raxml.mlTrees",
    output:
        outfile     = f"{full_dir}filtered_trees",
    script:
        "scripts/filter_iqtrees.py"

rule collect_best_llh:
    """
    This rule collects the likelihood values of the best tree as given at the end of the raxml-ng treesearch 
    and the likelihood of the same tree as evaluated by iqtree. 
    """
    input: 
        raxml_treesearch_log    = f"{full_file_path_raxml}.raxml.treesearch.log",
        iqtree_test_log         = f"{full_file_path_iqtree}.iqtree_tests.log",
    output:
        outfile     = f"{full_dir}llh_bestTree",
    run:
        # read final llh as computed in raxml-ng run
        with open(input.raxml_treesearch_log) as f:
            content = f.readlines()
        if not any("Final LogLikelihood" in l for l in content):
            raise ValueError(f"The given input file {input.raxml_treesearch_log} does not contain a final llh. Make sure the files are correct.") 
        for line in content:
            if "Final LogLikelihood" in line:
                _, llh = line.split(":")
                llh_raxml = llh.strip()
        
        # read final llh as reevaluated by iqtree
        with open(input.iqtree_test_log) as f:
            content = f.readlines()
        if not any("BEST SCORE FOUND" in l for l in content):
            raise ValueError(f"The given input file {input.iqtree_test_log} does not contain a best score. Make sure the files are correct.") 
        for line in content:
            if "BEST SCORE FOUND" in line:
                _, llh = line.split(":")
                llh_iqtree = llh.strip()
        
        # write the llhs to the output file
        # save also parameter settings in the same file
        with open(output.outfile, "w") as w:
            w.write(f"blmin:{wildcards.blmin}\nblmax:{wildcards.blmax}\n")
            w.write(f"Raxml-ng llh:{llh_raxml}\nIqtree llh:{llh_iqtree}\n")


