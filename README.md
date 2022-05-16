# Snakemake Pipeline

1. Clone this repo and cd into the snakemake directory
2. Create a conda environment using the given environment file by executing
    ```
    conda env create  --file environment.yml
    ```
3. Activate the conda environment

4. In the config.yaml file: 
   * specifiy the path to the input data and the output path by changing the variables `data.input` and `data.outdir`
   * Set the `parameters.model` to the substitution model you want to use. Refer to the manuals of RAxML-NG, IQ-Tree, and FastTree to set the correct string. 
   * Change the software execution paths for RAxML-NG, IQ-Tree, FastTree

5. In your terminal: check the execution graph and the commands snakemake will execute by running.
    ```
    snakemake -np
    ```

6. Finally, execute the snakemake file with 
    ```
    snakemake --cores [however many cores you have available]
    ```


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


# Branches
This repository consists of 4 branches:
1. master: This branch contains the pipeline we used to run RAxML-NG, IQ-Tree, and FastTree during Study 1 on Data Collection 1.
2. raxmlng_new_datasets: This branch contains the pipeline we used to verify our findings for the likelihood epsilon thresholds in RAxML-NG on Data Collection 2.
3. iqtree_new_datasets: This branch contains the pipeline we used to verify our findings for the likelihood epsilon threshold in IQ-Tree on Data Collection 2.
4. raxmlng_study2: This branch contains the pipeline we used to analyze the separated likelihood epsilon thresholds in RAxML-NG during Study 2 on Data Collection 2.