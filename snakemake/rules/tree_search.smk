rule raxml_pars_tree:
    input:
        msa             = config["data"]["input"],
    output:
        raxml_best_tree       = f"{full_file_path_raxml_pars}.raxml.bestTree",
        raxml_best_model      = f"{full_file_path_raxml_pars}.raxml.bestModel",
    params:
        # general params
        model           = config["parameters"]["model"]["raxml-ng"],
        # raxml-ng specific params
        threads         = config["parameters"]["raxml-ng"]["threads"],
        prefix          = full_file_path_raxml_pars,
    log:
        raxml_treesearch_log    = f"{full_file_path_raxml_pars}.raxml.treesearch.log",
    shell:
        "{raxml_command} " 
        "--msa {input.msa} "
        "--model {params.model} "
        "--prefix {params.prefix} "
        "--blmin {wildcards.blmin} "
        "--blmax {wildcards.blmax} "
        "--threads {params.threads} "
        "--seed {wildcards.seed} "
        "--tree pars{{1}} "
        "> {log} "
    
rule raxml_rand_tree:
    input:
        msa             = config["data"]["input"],
    output:
        raxml_best_tree       = f"{full_file_path_raxml_rand}.raxml.bestTree",
        raxml_best_model      = f"{full_file_path_raxml_rand}.raxml.bestModel",
    params:
        # general params
        model           = config["parameters"]["model"]["raxml-ng"],
        # raxml-ng specific params
        threads         = config["parameters"]["raxml-ng"]["threads"],
        prefix          = full_file_path_raxml_rand,
    log:
        raxml_treesearch_log    = f"{full_file_path_raxml_rand}.raxml.treesearch.log",
    shell:
        "{raxml_command} " 
        "--msa {input.msa} "
        "--model {params.model} "
        "--prefix {params.prefix} "
        "--blmin {wildcards.blmin} "
        "--blmax {wildcards.blmax} "
        "--threads {params.threads} "
        "--seed {wildcards.seed} "
        "--tree rand{{1}} "
        "> {log} "

rule collect_all_trees:
    input:
        pars_trees = expand(f"{full_file_path_raxml_pars}.raxml.bestTree", seed=pars_seeds, allow_missing=True),
        rand_trees = expand(f"{full_file_path_raxml_rand}.raxml.bestTree", seed=rand_seeds, allow_missing=True),
    output:
        all_trees = f"{full_file_path_raxml}.allTreesCollected"
    shell:
        "cat {input.pars_trees} {input.rand_trees} > {output.all_trees}"

rule collect_all_logs:
    input:
        pars_trees_logs = expand(f"{full_file_path_raxml_pars}.raxml.treesearch.log", seed=pars_seeds, allow_missing=True),
        rand_trees_logs = expand(f"{full_file_path_raxml_rand}.raxml.treesearch.log", seed=rand_seeds, allow_missing=True),
    output:
        all_logs = f"{full_file_path_raxml}.allTreesearchLogs"
    shell:
        "cat {input.pars_trees_logs} {input.rand_trees_logs} > {output.all_logs}"

rule save_best_tree_to_file:
    input:
        trees = rules.collect_all_trees.output.all_trees,
        logs = rules.collect_all_logs.output.all_logs,
    output:
        best_tree = f"{full_file_path_raxml}.bestTreeOfRun"
    script:
        "scripts/save_best_tree.py"

rule re_eval_best_tree:
    input:
        msa                 = config["data"]["input"],
        best_tree_of_run    = f"{full_file_path_raxml}.bestTreeOfRun"
    output:
        f"{full_file_path_raxml_eval}.raxml.log"
        f"{full_file_path_raxml_eval}.raxml.bestTree"
    params:
        model           = config["parameters"]["model"]["raxml-ng"],
        threads         = config["parameters"]["raxml-ng"]["threads"],
        prefix          = full_file_path_raxml_eval,
    log:
        f"{full_file_path_raxml_eval}.raxml.eval.log"
    shell:
        "{raxml_command} "
        "--eval "
        "--tree {input.best_tree_of_run} "
        "--msa {input.msa} "
        "--model {params.model} "
        "--prefix {params.prefix} "
        "--blmin {wildcards.blmin} "
        "--blmax {wildcards.blmax} "
        "--threads {params.threads} "
        "> {log} "

rule collect_all_eval_trees:
    input:
        eval_trees = expand(f"{full_file_path_raxml_eval}.raxml.bestTree", blmin_eval=blmin_opts, blmax_eval=blmax_opts, allow_missing=True)
    output:
        all_eval_trees = f"{full_file_path_raxml}.allEvalTreesCollected"
    shell:
        "cat {input.eval_trees} > {output.all_eval_trees}"

rule collect_all_eval_logs:
    input:
        eval_logs = expand(f"{full_file_path_raxml_eval}.raxml.eval.log", blmin_eval=blmin_opts, blmax_eval=blmax_opts, allow_missing=True)
    output:
        all_eval_logs = f"{full_file_path_raxml}.allEvalLogs"
    shell:
        "cat {input.eval_logs} > {output.all_eval_logs}"

rule save_best_eval_tree_to_file:
    input:
        trees = rules.collect_all_eval_trees.output.all_eval_trees,
        logs = rules.collect_all_eval_logs.output.all_eval_logs,
    output:
        best_tree = f"{full_file_path_raxml}.bestEvalTreeOfRun"
    script:
        "scripts/save_best_tree.py"