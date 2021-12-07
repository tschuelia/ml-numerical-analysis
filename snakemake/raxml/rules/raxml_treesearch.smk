rule raxml_pars_tree:
    input:
        msa             = config["data"]["input"],
    output:
        raxml_best_tree       = f"{full_file_path_raxml_pars}.raxml.bestTree",
        raxml_best_model      = f"{full_file_path_raxml_pars}.raxml.bestModel",
        raxml_start_tree      = f"{full_file_path_raxml_pars}.raxml.startTree",
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
        "--lh-epsilon-auto {wildcards.lheps_auto} "
        "--lh-epsilon-fast {wildcards.lheps_fast} "
        "--lh-epsilon-slow {wildcards.lheps_slow} "
        "--lh-epsilon-brlen-full {wildcards.lheps_full} "
        "--lh-epsilon-brlen-triple {wildcards.lheps_trip} "
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
        raxml_start_tree      = f"{full_file_path_raxml_rand}.raxml.startTree",
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
        "--lh-epsilon-auto {wildcards.lheps_auto} "
        "--lh-epsilon-fast {wildcards.lheps_fast} "
        "--lh-epsilon-slow {wildcards.lheps_slow} "
        "--lh-epsilon-brlen-full {wildcards.lheps_full} "
        "--lh-epsilon-brlen-triple {wildcards.lheps_trip} "
        "--threads {params.threads} "
        "--seed {wildcards.seed} "
        "--tree rand{{1}} "
        "> {output.raxml_log} "
