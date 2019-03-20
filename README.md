
# NET-seq and RNA-seq analysis pipeline

## description

An analysis pipeline for NET-seq or RNA-seq data with the following major steps:

- processing of the 5' molecular barcode, if present
- 3' adapter and quality trimming with [cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html)
    - if analyzing RNA-seq, also trim poly-A tails
- alignment with [Tophat2](https://ccb.jhu.edu/software/tophat/index.shtml)
- selection of unique mappers
- removal of PCR duplicates, if molecular barcodes were used
- summaries of quality statistics from [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
- summaries of library processing statistics
- generation of coverage tracks
- library size and spike-in normalization of coverage
- genome-wide scatterplots and correlations
- *ab initio* transcript annotation with [StringTie](https://ccb.jhu.edu/software/stringtie/)
- differential expression analysis over reference transcripts and discovered transcripts with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- classification of discovered transcripts into genomic categories
- data visualization (heatmaps and metagenes, with the option to separate data into clusters of similar signal)

## requirements

### required software

- Unix-like operating system (tested on CentOS 7.2.1511)
- Git
- [conda](https://conda.io/docs/user-guide/install/index.html)
- [build-annotations pipeline](https://github.com/winston-lab/build-annotations)

### required files

- FASTQ files of NET-seq libraries prepared as described in [our publication](https://doi.org/10.1016/j.molcel.2018.09.005) or RNA-seq libraries of the same library structure (legacy RNA-seq libraries sequenced from the 5'-end are also supported). Libraries with or without random hexamer molecular barcodes are both supported, with or without spikeins. The pipeline has only been tested using Illumina sequencing data. FASTQ files should be demultiplexed but not otherwise modified.

- FASTA files:
    - the 'experimental' genome
    - the spikein genome, if any samples have spikeins

- [BED6](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format annotation files:
    - ORF annotation
    - transcript annotation
    - optional: other annotations for data visualization (i.e. heatmaps and metagenes)

## instructions
**0**. If you haven't already done so, clone the separate ['build-annotations' pipeline](https://github.com/winston-lab/build-annotations), make a copy of the `config_template.yaml` file called `config.yaml`, and edit `config.yaml` as needed so that it points to the experimental genome FASTA file, ORF annotation BED file, and transcript annotation BED file to be used for the NET/RNA-seq pipeline. The 'build-annotations' pipeline will be used to create annotation files needed for classifying discovered transcripts into genomic categories.

**1**. Clone this repository.

```bash
git clone https://github.com/winston-lab/net-and-rna-seq.git
```

**2**. Create and activate the `snakemake_default` virtual environment for the pipeline using conda. The virtual environment creation can take a while. If you've already created the `snakemake_default` environment from another one of my pipelines, this is the same environment, so you can skip creating the environment and just activate it.

```bash
# navigate into the pipeline directory
cd net-and-rna-seq

# create the snakemake_default environment
conda env create -v -f envs/snakemake_default.yaml

# activate the environment
source activate snakemake_default

# to deactivate the environment
# source deactivate
```

**3**. Make a copy of the configuration file template `config_template.yaml` called `config.yaml`, and edit `config.yaml` to suit your needs.

```bash
# make a copy of the configuration template file
cp config_template.yaml config.yaml

# edit the configuration file
vim config.yaml    # or use your favorite editor
```

**4**. With the `snakemake_default` environment activated, do a dry run of the pipeline to see what files will be created.

```bash
snakemake -p --use-conda --dryrun
```

**5**. If running the pipeline on a local machine, you can run the pipeline using the above command, omitting the `--dryrun` flag. You can also use N cores by specifying the `--cores N` flag. The first time the pipeline is run, conda will create separate virtual environments for some of the jobs to operate in. Running the pipeline on a local machine can take a long time, especially for many samples, so it's recommended to use an HPC cluster if possible. On the HMS O2 cluster, which uses the SLURM job scheduler, entering `sbatch slurm_submit.sh` will submit the pipeline as a single job which spawns individual subjobs as necessary. This can be adapted to other job schedulers and clusters by adapting `slurm_submit.sh`, which submits the pipeline to the cluster, `slurm_status.sh`, which handles detection of job status on the cluster, and `cluster.yaml`, which specifies the resource requests for each type of job.

