rule re_eval_best_raxml_tree:
    input:
        msa                 = config["data"]["input"],
        best_tree_of_run    = f"{full_file_path_raxml}.bestTreeOfRun",
    output:
        log         = f"{full_file_path_raxml_eval}.raxml.log",
        best_tree   = f"{full_file_path_raxml_eval}.raxml.bestTree",
        eval_log    = f"{full_file_path_raxml_eval}.raxml.eval.log",
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
        "--lh-epsilon {wildcards.lh_eps_eval} "
        "--param-eps {wildcards.raxml_param_epsilon_eval} "
        "--brlen-smoothings {wildcards.raxml_brlen_smoothings_eval} "
        "--threads {params.threads} "
        "> {output.eval_log} "