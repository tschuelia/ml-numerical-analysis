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
        blmin_eval      = blmin_eval,
        blmax_eval      = blmax_eval,
        lh_eps_eval     = lh_eps_eval,
        model_param_epsilon_eval    = model_param_epsilon_eval,
        raxml_brlen_smoothings_eval = raxml_brlen_smoothings_eval,
        bfgs_fac_eval               = bfgs_fac_eval,
    log:
        f"{full_file_path_raxml_eval_pars}.snakelog"
    shell:
        "{raxml_command} "
        "--eval "
        "--tree {input.best_tree_of_run} "
        "--msa {input.msa} "
        "--model {params.model} "
        "--prefix {params.prefix} "
        "--blmin {params.blmin_eval} "
        "--blmax {params.blmax_eval} "
        "--lh-epsilon {params.lh_eps_eval} "
        "--param-eps {params.model_param_epsilon_eval} "
        "--brlen-smoothings {params.raxml_brlen_smoothings_eval} "
        "--bfgs-factor {params.bfgs_fac_eval} "
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
        model=config["parameters"]["model"]["raxml-ng"],
        threads=config["parameters"]["raxml-ng"]["threads"],
        prefix=full_file_path_raxml_eval_rand,
        blmin_eval=blmin_eval,
        blmax_eval=blmax_eval,
        lh_eps_eval=lh_eps_eval,
        model_param_epsilon_eval=model_param_epsilon_eval,
        raxml_brlen_smoothings_eval=raxml_brlen_smoothings_eval,
        bfgs_fac_eval=bfgs_fac_eval,
    log:
        f"{full_file_path_raxml_eval_rand}.snakelog"
    shell:
        "{raxml_command} "
        "--eval "
        "--tree {input.best_tree_of_run} "
        "--msa {input.msa} "
        "--model {params.model} "
        "--prefix {params.prefix} "
        "--blmin {params.blmin_eval} "
        "--blmax {params.blmax_eval} "
        "--lh-epsilon {params.lh_eps_eval} "
        "--param-eps {params.model_param_epsilon_eval} "
        "--brlen-smoothings {params.raxml_brlen_smoothings_eval} "
        "--bfgs-factor {params.bfgs_fac_eval} "
        "--threads {params.threads} "
        "> {output.eval_log} "