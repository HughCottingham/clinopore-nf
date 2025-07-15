# clinopore-nf

## Overview

An ONT-first bacterial isolate assembly pipeline using the best combination of **fully automated** tools according to the Trycycler(1) and Polypolish papers(2). Can be used with long reads alone or with optional short reads for polishing.

It takes a set of long reads and optional short reads and runs the following tools (all steps optional):

- filtlong (pre-assembly filtering) (https://github.com/rrwick/Filtlong)
- flye (assembly) (https://github.com/fenderglass/Flye)
- medaka (long-read polishing) (https://github.com/nanoporetech/medaka)
- polypolish (short-read polishing) (https://github.com/rrwick/Polypolish)
- polca (short-read polishing) (https://github.com/alekseyzimin/masurca)

## System Requirements

This pipeline was written for use in Linux operating systems and requires mamba/conda for installation. It has been tested on Ubuntu 18.04.5 and Rocky Linux 9.4 (Blue Onyx). Installation with mamba typically takes 4-5 minutes. Hybrid assembly of the demo readsets available at SRA accessions SRR13501341 (Illumina) and SRR34436406 (ONT) took 16 minutes on a Ubuntu 18.04.5 system using 8 CPUs and 16GB of RAM.

## Quickstart

```bash
#Open up a screen session

screen -R run_clinopore

#Clone the clinopore repository:

git clone https://github.com/HughCottingham/clinopore-nf.git && cd clinopore-nf

#Install the conda environments (skip this step if running on MASSIVE):

conda env create -f clinopore.yml -p /path/to/clinopore/conda/env
conda env create -f polca.yml -p /path/to/polca/conda/env

#Set run parameters:

nano nextflow.config

#Run the pipeline

conda activate /path/to/clinopore/conda/env
nextflow run clinopore.nf --clinopore_env /path/to/clinopore/conda/env --polca_env /path/to/polca/conda/env
```

## Usage

The most important part of running the pipeline is ensuring your read data is in the correct format and has its path set as an input parameter. All reads (short and long) should be specified in the `nextflow.config` file (`reads/*fastq.gz` by default). Long reads should be in the format `{isolate_id}.fastq.gz` while short reads should be in the format `{isolate_id}_1.fastq.gz`/`{isolate_id}_2.fastq.gz`.

There are several other options that can be set by modifying `nextflow.config`. Below is a summary of these options.

| Option                            | Description                                                       | Default           |
| ----                              | ----                                                              | ----              |
| `clinopore_env`           | clinopore conda environment (MUST set)     | '/path/to/clinopore/env'                 |
| `polca_env`           | polca conda environment (MUST set)     | '/path/to/polca/env'                 |
| `reads`                | Glob pattern matching input reads             | 'reads/*fastq.gz'
| `medaka_model`              | Medaka model (see medaka github documentation for details)   | 'r941_min_sup_g507'                |
| `filter_reads`                   | Optional read filtering step ('true' or 'false')     | 'true'                 |
| `long_read_polish`               | Optional long-read polishing step ('true' or 'false')     | 'true'                |
| `short_read_polish`               | Optional short-read polishing step ('true' or 'false')     | 'true'                |
| `outdir`                        | Output directory                                                   | 'assemblies'   |
| `read_qc`              | Optional read QC (Nanostats and fastqc) ('true' or 'false')     | 'true'                |
| `polishing_stats`              | Optional summary of bases added/subtracted at each polishing step ('true' or 'false')     | 'true'                |
| `failure_action`              | Option to ignore ('ignore') or terminate ('terminate') on errors (local execution only)     | 'terminate'                |
| `threads`           | Number of threads to use (local execution only)      | 16                |
| `max_retries`              | MASSIVE use only - How many times should slurm retry the job with higher CPUs, RAM and/or walltime?   | 3                |
| `queue_size`                    | MASSIVE use only - Maximum number of jobs to submit at once                         | 50                 |
| `slurm_account`              | MASSIVE use only (MUST set if running on MASSIVE)   | 'slurm_account'                |

## Output

If all steps of the pipeline are run, the final assembly files will be in the base of the output directory named `{isolate_id}_medaka_polypolish_polca.fasta`. If you don't run every step or are interested in the intermediates, the intermediate fasta files are in the directories `{output_dir}/flye`, `{output_dir}/medaka` and `{output_dir}/polypolish`. The assembly graph files produced by flye are also provided in the `{output_dir}/gfa` directory to allow for visualisation in a tool like Bandage. 

Contigs are sorted from largest to smallest, meaning the chromosome will typically be `contig_1`, the largest plasmid will be `contig_2` and so on. Fasta headers also include contig length in bases, depth and whether it is circular or not. 


## References

1. https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02483-z, https://github.com/rrwick/Trycycler
2. https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009802, https://github.com/rrwick/Polypolish

