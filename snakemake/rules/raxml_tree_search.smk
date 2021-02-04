rule raxml_pars_tree:
    input:
        msa             = config["data"]["input"],
    output:
        raxml_best_tree       = f"{full_file_path_raxml_pars}.raxml.bestTree",
        raxml_best_model      = f"{full_file_path_raxml_pars}.raxml.bestModel",
        raxml_log             = f"{full_file_path_raxml_pars}.raxml.treesearch.log",
    params:
        # general params
        model           = config["parameters"]["model"]["raxml-ng"],
        # raxml-ng specific params
        threads         = config["parameters"]["raxml-ng"]["threads"],
        prefix          = full_file_path_raxml_pars,
    log:
        f"{full_file_path_raxml_pars}.snakelog",
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
        "> {output.raxml_log} "
    
rule raxml_rand_tree:
    input:
        msa             = config["data"]["input"],
    output:
        raxml_best_tree       = f"{full_file_path_raxml_rand}.raxml.bestTree",
        raxml_best_model      = f"{full_file_path_raxml_rand}.raxml.bestModel",
        raxml_log             = f"{full_file_path_raxml_rand}.raxml.treesearch.log",
    params:
        # general params
        model           = config["parameters"]["model"]["raxml-ng"],
        # raxml-ng specific params
        threads         = config["parameters"]["raxml-ng"]["threads"],
        prefix          = full_file_path_raxml_rand,
    log:
        f"{full_file_path_raxml_rand}.snakelog",
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
        "> {output.raxml_log} "

rule collect_all_raxml_trees:
    input:
        pars_trees = expand(f"{full_file_path_raxml_pars}.raxml.bestTree", seed=pars_seeds, allow_missing=True),
        rand_trees = expand(f"{full_file_path_raxml_rand}.raxml.bestTree", seed=rand_seeds, allow_missing=True),
    output:
        all_raxml_trees = f"{full_file_path_raxml}.allTreesCollected"
    shell:
        "cat {input.pars_trees} {input.rand_trees} > {output.all_raxml_trees}"

rule collect_all_raxml_logs:
    input:
        pars_trees_logs = expand(f"{full_file_path_raxml_pars}.raxml.treesearch.log", seed=pars_seeds, allow_missing=True),
        rand_trees_logs = expand(f"{full_file_path_raxml_rand}.raxml.treesearch.log", seed=rand_seeds, allow_missing=True),
    output:
        all_raxml_logs = f"{full_file_path_raxml}.allTreesearchLogs"
    shell:
        "cat {input.pars_trees_logs} {input.rand_trees_logs} > {output.all_raxml_logs}"

rule save_best_raxml_tree_to_file:
    input:
        trees = rules.collect_all_raxml_trees.output.all_raxml_trees,
        logs = rules.collect_all_raxml_logs.output.all_raxml_logs,
    output:
        best_tree = f"{full_file_path_raxml}.bestTreeOfRun"
    params:
        program = "raxml"
    script:
        "scripts/save_best_tree.py"

rule re_eval_best_raxml_tree:
    input:
        msa                 = config["data"]["input"],
        best_tree_of_run    = f"{full_file_path_raxml}.bestTreeOfRun",
    output:
        f"{full_file_path_raxml_eval}.raxml.log",
        f"{full_file_path_raxml_eval}.raxml.bestTree",
        eval_log = f"{full_file_path_raxml_eval}.raxml.eval.log",
    params:
        model           = config["parameters"]["model"]["raxml-ng"],
        threads         = config["parameters"]["raxml-ng"]["threads"],
        prefix          = full_file_path_raxml_eval,
    log:
        f"{full_file_path_raxml_eval}.snakelog"
    shell:
        "{raxml_command} "
        "--eval "
        "--tree {input.best_tree_of_run} "
        "--msa {input.msa} "
        "--model {params.model} "
        "--prefix {params.prefix} "
        "--blmin {wildcards.blmin_eval} "
        "--blmax {wildcards.blmax_eval} "
        "--threads {params.threads} "
        "> {output.eval_log} "

rule collect_all_raxml_eval_trees:
    input:
        eval_trees = expand(f"{full_file_path_raxml_eval}.raxml.bestTree", blmin_eval=blmin_opts, blmax_eval=blmax_opts, allow_missing=True)
    output:
        all_eval_trees = f"{full_file_path_raxml}.allEvalTreesCollected"
    shell:
        "cat {input.eval_trees} > {output.all_eval_trees}"

rule collect_all_raxml_eval_logs:
    input:
        eval_logs = expand(f"{full_file_path_raxml_eval}.raxml.eval.log", blmin_eval=blmin_opts, blmax_eval=blmax_opts, allow_missing=True)
    output:
        all_eval_logs = f"{full_file_path_raxml}.allEvalLogs"
    shell:
        "cat {input.eval_logs} > {output.all_eval_logs}"

rule save_best_raxml_eval_tree_to_file:
    input:
        trees = rules.collect_all_raxml_eval_trees.output.all_eval_trees,
        logs = rules.collect_all_raxml_eval_logs.output.all_eval_logs,
    output:
        best_tree = f"{full_file_path_raxml}.bestEvalTreeOfRun"
    params:
        program = "raxml"
    script:
        "scripts/save_best_tree.py"