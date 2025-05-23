# RNA_seq analysis pipeline
The scripts and pipelines used to analyze data for () manuscript.

Docker images that are called in this script can be found here:
* [general_docker](https://hub.docker.com/repository/docker/kwellswrasman/rnaseq_general/general) (Version 1)
* [r_docker](https://hub.docker.com/repository/docker/kwellswrasman/rnaseq_r/) (Version 1)


## Running the pipeline

A helper script to create the sample sheet and to start creating the rmats and deseq sample info can be run with 
```bash
python src/scripts/make_sample_sheets.py
```

A few notes about this script
1. The default fastq directory is `raw_data`, however you can provide your own with `-d`.
2. The default output directory is `files`, however you can provide your own with `-o`.
3. The `samples.tsv` file will be fully ready to use with the snakemake pipeline.
4. The `sample_info.csv` file will be made with samples and columns added in, but you will need to fill it out with your own sample information.
5. The `rmats_sample_info.tsv` file will be made with samples and columns added in, but you will need to fill it out with your own sample information.
6. A warning is printed when running the script to remind you that these two files are not complete.

## Snakemake

A snakemake pipeline that can be used to run bulk RNA-seq analysis. Can chose between cutadapt, bbduk or no adapter trimming. Outputs fastqc summary files, star summary files, and a counts matrix, and an html file containing quality control images. All of the output files can be analyzed using the rmd script

Writen by Kristen Wells

To use:

1. Download and install miniconda3: For Linux
```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh bash Miniconda3-latest-Linux-x86_64.sh
```
2. Install Snakemake:
```bash
conda install snakemake -c bioconda -c conda-forge
```

This pipeline has been configured to run on a slurm scheduler. For more information on running snakemake on different cluster configerations see the [executor help page](https://snakemake.github.io/snakemake-plugin-catalog/index.html)

Changing the profile should be the only major update to move to a different scheduler, however, I have added information for custom log files to each rule that is specific to slurm.

3. Update the config file (config.yaml) 
>* SAMPLE_TABLE: path to a file consisting of at least two columns.
>   * The first column should be titled `sample` and contain the name of the sample.
>   * The second should be titled `fastq1` and contain the path to the associated fastq file
>   * The third *optional* column should be titled `fastq2` and be used if you have paired end data. This should contain the path to the associated read 2 fastq file.
>   * The fourth *optional* column is named `spike-in` with TRUE/FALSE values per sample indicating if spikeins were used. If this column is not included it will be assumed that spike-ins were not included.
>* PROJECT: The name of the project, will be the name of the output counts matrix
>* GENOME: The path to the genome directory created by star
>* COMBINED_GENOME: *Optional*, only include this if spike-ins were used. This is the genome of both the main organism and spike-in organism
>* GTF: The path to the GTF associated with the star genome
>* COMBINED_GTF: *Optional*, only include with spike-ins, the path to the combined GTF file
>* RESULTS: The path to the results directory
>* ADAPTORS: *Optional*, only include if using bbduk for adaptor trimming. Path to the adaptors file in the bbtools package
>* TRIM_METHOD: What trim method to use. Can be "bbduk", "cutadapt", or "no"
>* PE and SE: extra paramamters for all of the jobs run

4. Update snakecharmer.sh to your specific cluster specs. 

5. submit the job using `sbatch snakecharmer.sh`

## R analysis
The file `src/scripts/DESeq_analysis.Rmd` was used for RNA-seq analysis. Figures were made with `src/scripts/figures.R`

