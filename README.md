# clinopore-nf

An ONT-first assembly pipeline using the best combination of fully automated tools according to the Trycycler(1) and Polypolish papers(2). Can be used with long reads alone or with optional short reads for polishing. If completing the manual steps required for running trycycler is not feasible, this is the next best option in many cases. 

It takes a set of long reads and optional short reads and runs the following tools:

- filtlong (pre-assembly filtering)
- flye (assembly)
- medaka (long-read polishing)
- polypolish (short-read polishing)
- polca (short-read polishing)

## Quickstart

```bash
#Open up a screen session

screen -R run_clinopore

#Clone the clinopore repository:

git clone https://github.com/HughCottingham/clinopore-nf.github&& cd clinopore-nf

#Install the conda environments (skip this step if running on MASSIVE):

conda create -f clinopore.yml
conda create -f polca.nf
conda env list # get the paths of installed envs - will need in run step

#Set input parameters:

nano nextflow.config

#Run the pipeline

nextflow run clinopore.nf --clinopore_env /path/to/clinopore/conda/env --polca_env /path/to/polca/conda/env
```

## Usage

The most important part of running the pipeline is ensuring your read data is in the correct format and has its path set as an input parameter. All reads (short and long) should be in the same specified directory (`reads/*fastq.gz` by default). This can be set by editing `nextflow.config`. Long reads should be in the format `{isolate_id}.fastq.gz` while short reads should be in the format `{isolate_id}_1.fastq.gz`/`{isolate_id}_2.fastq.gz`.

There are several other options that can be set by modifying `nextflow.config`. Below is a summary of these options.

| Option                            | Description                                                       | Default           |
| ----                              | ----                                                              | ----              |
| `reads`                | Glob pattern matching input reads               | 200               | 'reads/*fastq.gz'
| `outdir`                        | Output directory                                                   | 'assemblies'   |
| `run_filtlong`                   | Optional read filtering step ('true' or 'false')     | 'true'                 |
| `run_medaka`               | Optional long-read polishing step ('true' or 'false')     | 'true'                |
| `run_polypolish`               | Optional short-read polishing step ('true' or 'false')     | 'true'                |
| `run_polca`              |Optional short-read polishing step ('true' or 'false')     | 'true'               |
| `polca_env`                    | polca conda environment (MUST set)     | '/path/to/polca/env'                | 
| `clinopore_env`           | clinopore conda environment (MUST set)     | '/path/to/clinopore/env'                 |
| `max_retries`              | MASSIVE use only - How many times should slurm retry the job with higher CPUs, RAM and/or walltime?   | 3                |
| `queue_size`                    | MASSIVE use only - Maximum number of jobs to submit at once                         | 1000                 |
| `processors`           | Number of processors to use      | 4                |
| `slurm_account`              | MASSIVE use only (MUST set if running on MASSIVE)   | 'slurm_account'                |

## Output

If all steps of the pipeline are run, the final assembly files will be in the base of the output directory named `{isolate_id}_medaka_polypolish_polca.fasta`. If you don't run every step or are interested in the intermediates, the intermediate fasta files are in the directories `{output_dir}/flye`, `{output_dir}/medaka` and `{output_dir}/polypolish`. The assembly graph files produced by flye are also provided in the `{output_dir}/gfa` directory to allow for visualisation in a tool like Bandage. 

Contigs are sorted from largest to smallest, meaning the chromosome will typically be `contig_1`, the largest plasmid will be `contig_2` and so on. Fasta headers also include contig length in bases, depth and whether it is circular or not. 

## References

1. https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02483-z, https://github.com/rrwick/Trycycler
2. https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009802, https://github.com/rrwick/Polypolish
