rule reevaluate_raxml_pars_tree:
    input:
        msa                 = config["data"]["input"],
        best_tree_of_run    = f"{full_file_path_raxml_pars}.raxml.bestTree",
    output:
        log         = f"{full_file_path_raxml_eval_pars}.raxml.log",
        best_tree   = f"{full_file_path_raxml_eval_pars}.raxml.bestTree",
        eval_log    = f"{full_file_path_raxml_eval_pars}.raxml.eval.log",
    params:
        model           = config["parameters"]["model"]["raxml-ng"],
        threads         = config["parameters"]["raxml-ng"]["threads"],
        prefix          = full_file_path_raxml_eval_pars,
    log:
        f"{full_file_path_raxml_eval_pars}.snakelog"
    shell:
        "{raxml_command} "
        "--eval "
        "--tree {input.best_tree_of_run} "
        "--msa {input.msa} "
        "--model {params.model} "
        "--prefix {params.prefix} "
        "--lh-epsilon-auto 0.1 "
        "--lh-epsilon-fast 0.1 "
        "--lh-epsilon-slow 0.1 "
        "--lh-epsilon-brlen-full 0.1 "
        "--lh-epsilon-brlen-triple 0.1 "
        "--threads {params.threads} "
        "> {output.eval_log} "


rule reevaluate_raxml_rand_tree:
    input:
        msa                 = config["data"]["input"],
        best_tree_of_run    = f"{full_file_path_raxml_rand}.raxml.bestTree",
    output:
        log         = f"{full_file_path_raxml_eval_rand}.raxml.log",
        best_tree   = f"{full_file_path_raxml_eval_rand}.raxml.bestTree",
        eval_log    = f"{full_file_path_raxml_eval_rand}.raxml.eval.log",
    params:
        model   = config["parameters"]["model"]["raxml-ng"],
        threads = config["parameters"]["raxml-ng"]["threads"],
        prefix  = full_file_path_raxml_eval_rand,
    log:
        f"{full_file_path_raxml_eval_rand}.snakelog"
    shell:
        "{raxml_command} "
        "--eval "
        "--tree {input.best_tree_of_run} "
        "--msa {input.msa} "
        "--model {params.model} "
        "--prefix {params.prefix} "
        "--lh-epsilon-auto 0.1 "
        "--lh-epsilon-fast 0.1 "
        "--lh-epsilon-slow 0.1 "
        "--lh-epsilon-brlen-full 0.1 "
        "--lh-epsilon-brlen-triple 0.1 "
        "--threads {params.threads} "
        "> {output.eval_log} "