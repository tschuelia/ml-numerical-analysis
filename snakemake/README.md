# Snakemake Pipeline

1. Create a conda environment using the given environment file by executing
    ```
    conda env create --name bioinfo --file environment.yml
    ```
2. Activate the conda environment

3. Download the data you want to use, for example the [raxml-ng tutorial dataset](https://github.com/amkozlov/ng-tutorial) into a folder named ```data```

4. In the config.yaml file: specifiy the path to the input data and the output path by changing the variables ```data.input``` and ```data.outdir``` . 

5. In your terminal: check the execution graph and the commands snakemake will execute by running.
    ```
    snakemake -np
    ```

6. Finally execute the snakemake file with 
    ```
    snakemake --cores [however many cores you have available]
    ```

## Parameter Changes
In case you want to add new variable parameters to a run there are a few things that need updating:
* ```database.py```: add the new parameters to the Run definition
* ```save_results_to_database.py::create_Run```: update the implementation to also add the new params to the run object

<hr>


# Working on a Cluster: 
## Job submission
There are two ways to submit the snakemake job on the cluster.
The first is using a sbatch file and the second is using the power of snakemake. Using sbatch is easier to setup but using snakemake is more powerful since tasks are distributed more fine grained (at least I think so :-)).
### Slurm Jobs
Create a sbatch file ```job.sh```:

    ```
    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=2
    #SBATCH -t 00:20:00
    #SBATCH -p normal
    #SBATCH --mem=60000
    #SBATCH --chdir=aa_rokasA1_10@4
    #SBATCH -o out.txt
    #SBATCH -e err.txt

    *normal bash here*
    snakemake --cores all
    ```
and submit it using ```sbatch job.sh```


### Using snakemake:
Snakemake uses profiles to run on different machines. To use slurm do the following setup (according to [these instructions](http://bluegenes.github.io/Using-Snakemake_Profiles/)):

1. ssh to the cluster and run 
    ``` 
    mkdir -p ~/.config/snakemake && cd ~/.config/snakemake 
    ```
2. In this dir run 
    ```
        cookiecutter https://github.com/Snakemake-Profiles/slurm.git
    ```
    Leave all setup prompts empty for now (just press enter).
3. In ```~/.config/snakemake``` should now be ```slurm``` folder containing
    ```
    config.yaml
    slurm-jobscript.sh
    slurm-status.py
    slurm-submit.py
    slurm_utils.py  
    ```
4. Open ```config.yaml``` and add the snakemake settings you wish to use. Mine looks like this:
    ```
    restart-times: 3
    jobscript: "slurm-jobscript.sh"
    cluster: "slurm-submit.py"
    cluster-status: "slurm-status.py"
    max-jobs-per-second: 1
    max-status-checks-per-second: 10
    local-cores: 1
    latency-wait: 60
    use-conda: True
    jobs: 10
    rerun-incomplete: True
    printshellcmds: True
    ```
5. Now we have to define the resource utilization. For this create a file in the same dir called ```cluster_config.yml``` and add the desired default resources. Mine looks like this:
    ```
    __default__:
        partition: partition # the partition you use
        mail-user: your@email.com 
        time: 60 # default time in minutes
        nodes: 1 
        cpus-per-task: 1 
        threads-per-core: 1
        output: out.txt
        error: err.txt
        mem: 1GB # default memory
    ```
If your rules do not specify custom ```resources``` snakemake will use these settings.

6. Open ```slurm-submit.py``` and edit the ```CLUSTER_CONFIG``` constant so it points to your ```cluster_config.yml``` file.

7. Open the folder containing the ```Snakefile``` and run 
    ```
    snakemake --profile slurm
    ```

**Note:** The above settings are use relatively little resources. Some of the rules, for example ```rules.tree_search.smk:raxml_tree``` need more resources. You can override the default settings by specifying resources in each rule. For some rules these options are already set. If you change the dataset, these options must probably be varied.



## Useful stuff
### Open ssh tmux session in iTerm 
 ssh -tt forhlr2 "tmux -CC new -A -s main"

 raxml-ng --parse --msa data/covid_data/fmsan/data/covid_edited.fasta --model GTR+G --prefix test

### Show Cluster Usage
```
sacct --format="JobID,JobName,NTasks,NCPUS,NNodes,Timelimit,AllocCPUS,State,Submit,End" 
```
