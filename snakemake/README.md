# Snakemake Pipeline

1. Create a conda environment using the given environment file by executing
    ```
    conda env create --name bioinfo --file environment.yml
    ```
2. Activate the conda environment

3. Download the IQ-Tree binaries into a folder ```software/iqtree```

4. Download the data you want to use, for example the [raxml-ng tutorial dataset](https://github.com/amkozlov/ng-tutorial) into a folder named ```data```

5. In the config.yaml file: specifiy the path to the input data and the output path by changing the variables ```data.input``` and ```data.outdir``` . 

6. In your terminal: check the execution graph and the commands snakemake will execute by running.
    ```
    snakemake -np
    ```

7. Finally execute the snakemake file with 
    ```
    snakemake --cores [however many cores you have available]
    ```

