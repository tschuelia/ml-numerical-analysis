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
        "--lh-epsilon {wildcards.lh_eps} "
        "--param-eps {wildcards.model_param_epsilon} "
        "--brlen-smoothings {wildcards.raxml_brlen_smoothings} "
        "--spr-lheps {wildcards.spr_lh_eps} "
        "--bfgs-factor {wildcards.bfgs_fac}"
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
        "--lh-epsilon {wildcards.lh_eps} "
        "--param-eps {wildcards.model_param_epsilon} "
        "--brlen-smoothings {wildcards.raxml_brlen_smoothings} "
        "--spr-lheps {wildcards.spr_lh_eps} "
        "--bfgs-factor {wildcards.bfgs_fac}"
        "--seed {wildcards.seed} "
        "--tree rand{{1}} "
        "> {output.raxml_log} "
