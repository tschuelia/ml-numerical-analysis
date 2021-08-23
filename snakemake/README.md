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
    #SBATCH -o out.txt
    #SBATCH -e err.txt

    *normal bash here*
    snakemake --cores all
    ```
and submit it using ```sbatch job.sh```


### Using snakemake:
Follow the instructions [here](https://github.com/tschuelia/snakemake-on-slurm-clusters).

## Useful stuff
### Open ssh tmux session in iTerm 
 ssh -tt forhlr2 "tmux -CC new -A -s main"

 raxml-ng --parse --msa data/covid_data/fmsan/data/covid_edited.fasta --model GTR+G --prefix test

### Show Cluster Usage
```
sacct --format="JobID,JobName,NTasks,NCPUS,NNodes,Timelimit,AllocCPUS,State,Submit,End" 
```
