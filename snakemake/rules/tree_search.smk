rule raxml_tree:
    input:
        msa             = config["data"]["input"],
    output:
        best_tree       = f"{full_file_path}.raxml.bestTree",
        best_model      = f"{full_file_path}.raxml.bestModel",
        ml_trees        = f"{full_file_path}.raxml.mlTrees",
    params:
        threads         = config["parameters"]["raxml-ng"]["threads"],
        seed            = config["parameters"]["raxml-ng"]["seed"],
        model           = config["parameters"]["raxml-ng"]["model"],
        num_pars_trees  = config["parameters"]["raxml-ng"]["num_pars_trees"],
        num_rand_trees  = config["parameters"]["raxml-ng"]["num_rand_trees"],
        prefix          = full_file_path,
    log:
        f"{full_file_path}.raxml.log",
    shell:
        "{raxml_command} --msa {input} --model {params.model} --prefix {params.prefix} "
        "--blmin {wildcards.blmin} --blmax {wildcards.blmax} "
        "--threads {params.threads} --seed {params.seed} "
        "--tree pars{{{params.num_pars_trees}}},rand{{{params.num_rand_trees}}} "
        "> {log} "