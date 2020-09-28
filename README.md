# CGI\_de\_novo\_real\_Pipeline

This pipeline is for Sentieon's assembly software for use with real sample data.

## Directory Contents

- Snakefile
    - This snakefile contains all the rules for running Sentieon's de novo assembly software using their seed library pipeline and real data
    - This includes
        - de novo assembly
        - variant calling based on different ploidies
        - quast analyses against a reference
        - phase block N50 calculations
        - Scaffolding of contigs
- config.yaml
    - The config file contains all the parameters necessary for modular assembly and evaluation
    - More info can be found in the comments and below
- phase\_block\_n50.py
    - Script to calculate the phase block N50
    - (1) Sum the size of all contigs by PID for all PID >= 0. 
    - (2) For each contig with a PID equal to -1, add the length of this contig to the larger adjacent PID within the scaffold.
    - (3) Order the PIDs from smallest to largest and calculate the phase N50 as usual.)))
- scaffold\_contigs.py
    - Script to scaffold contigs
    - This combines contigs with the same scaffold id (SID)
    - It inserts Ns such that they equal the GAP field of the contig

## Running the pipeline

Link the fastq files you'll be working with as `split_read.N.fq.gz` wherever you want to run the assembly.
I've run these typically in the datas directory.
An example would be `/research/rv-02/Projects/Project_stLFR/BGI/pipeline_analysis_for_combined_lanes/ZY20180118lib-O-NA1287850M50p11C-16_combine/Sentieon_CG_202003.18_Assembly/`

Then copy `config.yaml` to the same directory and update as appropriate.
I have snakemake installed as a virtual env.
You can source it using `source /home/eanderson/Virtual_Envs/SnakeMake/bin/activate`
From their run some version of the following.

```
# -j specifies the number of threads to use
# -k specifies keep going with independent jobs after an error
# -s specifies the location of the snakefile
# you can also copy the snakefile to the current directory and omit -s
snakemake -j 20 -k -s /research/rv-02/home/eanderson/CGI_de_novo_real_Pipeline/Snakefile 2>&1 | tee snakemake.err.txt
```

## Modifying config.yaml

- samples
    - fastq: `[split_read.1.fq.gz, split_read.2.fq.gz]`
        - modify this with the path to the split fastqs you want to use for Assembly
    - date:
        - I change this to the date run when I remember, but it doesn't affect the pipeline
- software
    - These shouldn't have to be changed too much
- params
    - sentieon_install: `path/to/newest/sentieon/software`
        - This should be updated when new software becomes available
    - seed_lib: `path/to/seed/lib`
        - This should be updated when new seed libs are developed
    - read_len: `200`
        - change as appropriate for the sample
    - is_female: `True`
        - sample dependent
- benchmark
    - benchmark\_snp: `/path/to/the/appropriate/snp/callset.vcf.gz`
        - should only have to be modified if you're benchmarking a sample that isn't HG001
    - benchmark\_indel: `/path/to/the/appropriate/indel/callset.vcf.gz`
    - bedfile: `/path/to/high/confidence/bedfile.bed`
    - all of the above are organized under `/research/rv-02/home/eanderson/Resources_And_DBs/hg38_GIAB/HG00N/<files>`
    - ref\_sdf: `/path/to/ref.sdf`
        - Should only need to be changed when the reference changes and you're doing benchmarking
- threads
    - These are various tools used by the pipeline
    - Specify the threads as you feel is appropriate
    - If you give the main snakemake command (run\_snakemake.sh) fewer threads, all tools will be capped at that number
