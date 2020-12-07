# Snakemake Playground

1. Create a conda environment using the given environment file by executing
    ```
    conda env create --name bioinfo --file environment.yml
    ```
2. Activate the conda environment

3. Download the [raxml-ng tutorial dataset](https://github.com/amkozlov/ng-tutorial) into a folder named ```data```

4. Check the execution graph by running 
    ```
    snakemake -n
    ```

5. Execute the snakemake file with 
    ```
    snakemake --cores 1
    ```

