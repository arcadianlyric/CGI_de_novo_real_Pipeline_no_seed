# Snakemake Pipeline for real de novo datasets
# Start from absolute scratch
# So need to make the assembly using the newest software (config file)
# prepend commands with export Sentieon stuff.

configfile: "config.yaml"
# These export the sentieon environment variables for all commands
shell.prefix('export SENTIEON_INSTALL={}; export SENTIEON_LICENSE={}; '.format(config['params']['sentieon_install'], config['params']['sentieon_license']))

# Rule for run all, could be shortened with a function
# These are the targets we want to make
rule run_all:
    input:
        "Quast_Contigs/quast.log",
        "hap0_phase_n50.txt",
        "hap1_phase_n50.txt",
        "Assembly/LinkedReads.contig_0.fasta",
        "Assembly/LinkedReads.contig_1.fasta",
        "Diploid_Vars/minimap_contigs0.bam",
        "Diploid_Vars/minimap_contigs1.bam",
        "Diploid_Vars/minimap_contigs0.bed",
        "Diploid_Vars/minimap_contigs1.bed",
        "Diploid_Vars/minimap_contig_phased_snp.vcf.gz",
        "Diploid_Vars/minimap_contig_phased_indel.vcf.gz",
        "Diploid_Vars/minimap_contig_diploid.bed",
        "Diploid_Vars/minimap_contig_haploid.bed",
        "Diploid_Vars/minimap_contig_missing.bed",
        "Diploid_Vars/minimap_contig_snp_haploid_vars/summary.txt",
        "Diploid_Vars/minimap_contig_snp_diploid_vars/summary.txt",
        "Diploid_Vars/minimap_contig_snp_missing_vars/summary.txt",
        "Diploid_Vars/minimap_contig_snp_missing_vars/summary.txt",
        "Diploid_Vars/minimap_contig_indel_haploid_vars/summary.txt",
        "Diploid_Vars/minimap_contig_indel_diploid_vars/summary.txt",
        "Diploid_Vars/minimap_contig_indel_missing_vars/summary.txt",
        "Diploid_Vars/minimap_squashed_contig_snp_haploid_vars/summary.txt",
        "Diploid_Vars/minimap_squashed_contig_snp_diploid_vars/summary.txt",
        "Diploid_Vars/minimap_squashed_contig_snp_missing_vars/summary.txt",
        "Diploid_Vars/minimap_squashed_contig_indel_haploid_vars/summary.txt",
        "Diploid_Vars/minimap_squashed_contig_indel_diploid_vars/summary.txt",
        "Diploid_Vars/minimap_squashed_contig_indel_missing_vars/summary.txt",
        "Scaffolds/Quast_Scaffolds/quast.log",
        "Scaffolds/Diploid_Vars/minimap_scaffolds0.bam",
        "Scaffolds/Diploid_Vars/minimap_scaffolds1.bam",
        "Scaffolds/Diploid_Vars/minimap_scaffolds0.bed",
        "Scaffolds/Diploid_Vars/minimap_scaffolds1.bed",
        "Scaffolds/Diploid_Vars/minimap_scaffold_phased_snp.vcf.gz",
        "Scaffolds/Diploid_Vars/minimap_scaffold_phased_indel.vcf.gz",
        "Scaffolds/Diploid_Vars/minimap_scaffold_snp_haploid_vars/summary.txt",
        "Scaffolds/Diploid_Vars/minimap_scaffold_snp_diploid_vars/summary.txt",
        "Scaffolds/Diploid_Vars/minimap_scaffold_snp_missing_vars/summary.txt",
        "Scaffolds/Diploid_Vars/minimap_scaffold_indel_haploid_vars/summary.txt",
        "Scaffolds/Diploid_Vars/minimap_scaffold_indel_diploid_vars/summary.txt",
        "Scaffolds/Diploid_Vars/minimap_scaffold_indel_missing_vars/summary.txt",
        "Scaffolds/Diploid_Vars/minimap_squashed_scaffold_snp_haploid_vars/summary.txt",
        "Scaffolds/Diploid_Vars/minimap_squashed_scaffold_snp_diploid_vars/summary.txt",
        "Scaffolds/Diploid_Vars/minimap_squashed_scaffold_snp_missing_vars/summary.txt",
        "Scaffolds/Diploid_Vars/minimap_squashed_scaffold_indel_haploid_vars/summary.txt",
        "Scaffolds/Diploid_Vars/minimap_squashed_scaffold_indel_diploid_vars/summary.txt",
        "Scaffolds/Diploid_Vars/minimap_squashed_scaffold_indel_missing_vars/summary.txt"


# Rule for read clustering
# mapping reads ot the seed library
rule read_clustering:
    input:
        expand("{sample}", sample=config['samples']['fastq'])
    output:
        "Assembly/cg50X_SentieonIter0.bam"
    log:
        "Assembly/run.log"
    threads:
        config['threads']['assembly']
    params:
        seed_lib = config['params']['seed_lib'],
        readgroup = r'@RG\tID:{0}\tSM:{0}\tPL:{1}'.format(config['samples']['date'],
                                                          config['params']['platform'])
    shell:
        "($SENTIEON_INSTALL/bin/sentieon bwa mem -R '{params.readgroup}' "
            "-M -t {threads} -K 10000000 -Y {params.seed_lib} {input} || echo -n 'BWA error') | "
            "extract | "
            "$SENTIEON_INSTALL/bin/sentieon util sort -t {threads} -o {output} --sam2bam - &> {log}"


# Sentieon's software uses a pcr indel model to limit false positives
rule indel_modeling:
    input:
        "Assembly/cg50X_SentieonIter0.bam"
    output:
        "Assembly/LinkedReads.pcr.model"
    log:
        "Assembly/LinkedReadsRepeatStat.log"
    threads:
        config['threads']['assembly']
    params:
        seed_lib = config['params']['seed_lib']
    shell:
        "$SENTIEON_INSTALL/bin/sentieon driver -t {threads} "
            "-r {params.seed_lib} -i {input} "
            "--algo RepeatStat {output} &> {log}"


# First assembly step
rule linked_reads_assemble:
    input:
        bam = "Assembly/cg50X_SentieonIter0.bam",
        model = "Assembly/LinkedReads.pcr.model"
    output:
        bam = "Assembly/LinkedReads.hap.bam"
    log:
        "Assembly/LinkedReadsAssemble.log"
    threads:
        config['threads']['assembly']
    params:
        seed_lib = config['params']['seed_lib'],
        read_len = config['params']['read_len']
    shell:
        "$SENTIEON_INSTALL/bin/sentieon driver -t {threads} "
            "-r {params.seed_lib} -i {input.bam} "
            "--algo LinkedReadsAssemble --assemble_read_len {params.read_len} "
            "--pcr_indel_model {input.model} "
            "--coalign_hap 1 --min_map_qual 1 "
            "{output.bam} &> {log}"


# First graph solving step
rule linked_reads_solve:
    input:
        clusters = "Assembly/cg50X_SentieonIter0.bam",
        bam = "Assembly/LinkedReads.hap.bam"
    output:
        graph = "Assembly/LinkedReads.hap.graph",
        reads = "Assembly/LinkedReads.hap.reads",
        filtered = "Assembly/LinkedReads.filtered.reads",
        rescue = "Assembly/LinkedReads.rescue.data"
    log:
        "Assembly/LinkedReadsSolve.log"
    threads:
        config['threads']['assembly']
    params:
        seed_lib = config['params']['seed_lib'],
        range_bc = config['params']['range_barcode']
    shell:
        "$SENTIEON_INSTALL/bin/sentieon driver -t {threads} "
            "-r {params.seed_lib} -i {input.clusters} "
            "--algo LinkedReadsSolve --hap_bam {input.bam} "
            "--hap_graph {output.graph} --reads_in_hap {output.reads} "
            "--reads_filtered {output.filtered} "
            "--rescue_data {output.rescue} "
            "--range_barcode {params.range_bc} &> {log}"


# rescue step for unassembled reads and contigs
rule linked_reads_rescue:
    input:
        bam = "Assembly/cg50X_SentieonIter0.bam",
        rescue = "Assembly/LinkedReads.rescue.data",
        reads = "Assembly/LinkedReads.hap.reads",
        filtered = "Assembly/LinkedReads.filtered.reads"
    output:
        "Assembly/LinkedReads.rescue.bam"
    log:
        "Assembly/LinkedReadsRescue.log"
    threads:
        config['threads']['assembly']
    shell:
        "$SENTIEON_INSTALL/bin/sentieon driver -t {threads} "
            "-i {input.bam} --algo LinkedReadsRescue "
            "--bam_compression 1 --rescue_data {input.rescue} "
            "--reads_in_hap {input.reads} "
            "--reads_filtered {input.filtered} "
            "{output} &> {log}"


# Second assembly step/scaffolding
rule linked_reads_bridge:
    input:
        bam = "Assembly/LinkedReads.rescue.bam",
        data = "Assembly/LinkedReads.rescue.data"
    output:
        "Assembly/LinkedReads.bridge.ret"
    log:
        "Assembly/LinkedReadsBridge.log"
    threads:
        config['threads']['assembly']
    params:
        seed_lib = config['params']['seed_lib']
    shell:
        "$SENTIEON_INSTALL/bin/sentieon driver -t {threads} "
            "-r {params.seed_lib} -i {input.bam} "
            "--algo LinkedReadsBridge --rescue_data {input.data} "
            "{output} &> {log}"


# Second graph solving step
rule linked_reads_solve_2:
    input:
        bam = "Assembly/cg50X_SentieonIter0.bam",
        hap = "Assembly/LinkedReads.hap.bam",
        graph = "Assembly/LinkedReads.hap.graph",
        bridge = "Assembly/LinkedReads.bridge.ret"
    output:
        contigs = expand("Assembly/LinkedReads.contig_{hap}.fasta", hap=[0,1]),
        phase_contigs = expand("Assembly/LinkedReads.contig.phase_{hap}.fasta", hap=[0,1]),
        barcodes = "Assembly/LinkedReads.phase.barcodes"
    log:
        "Assembly/LinkedReadsSolve1.log"
    threads:
        config['threads']['assembly']
    params:
        jaccard = config['params']['min_jaccard'],
        seed_lib = config['params']['seed_lib'],
        range_bc = config['params']['range_barcode']
    shell:
        "$SENTIEON_INSTALL/bin/sentieon driver -t {threads} "
            "-r {params.seed_lib} -i {input.bam} "
            "--algo LinkedReadsSolve --hap_bam {input.hap} "
            "--hap_graph {input.graph} --bridge_ret {input.bridge} "
            "--contig_output Assembly/LinkedReads.contig "
            "--min_scaffold_jaccard {params.jaccard} "
            "--range_barcode {params.range_bc} "
            "--phase_barcodes {output.barcodes} &> {log}"


# run quast on contig assemblies
rule run_quast_contigs:
    input:
        expand("Assembly/LinkedReads.contig_{hap}.fasta", hap=[0,1])
    output:
        "Quast_Contigs/quast.log"
    threads:
        config['threads']['quast']
    params:
        quast = config['software']['quast'],
        min_contig = config['params']['min_contig'],
        ref = config['params']['ref']
    shell:
        "{params.quast} {input} --eukaryote "
            "-R {params.ref} "
            "--threads {threads} "
            "--no-snps "
            "--fragmented "
            "--min-contig {params.min_contig} "
            "-o $(dirname {output}) "
            "-s "


# get contig N50 using python script
rule contig_phase_n50:
    input:
        "Assembly/LinkedReads.contig.phase_{hap}.fasta"
    output:
        "hap{hap}_phase_n50.txt"
    params:
        phase_script = config['software']['phase_script']
    shell:
        "{params.phase_script} {input} > {output}"


# Generate paf files using mimimap
rule run_minimap_paf:
    input:
        fasta = "Assembly/LinkedReads.contig_{hap}.fasta",
        ref = config['params']['ref']
    output:
        "Diploid_Vars/minimap_contigs{hap}.paf.gz"
    threads:
        config['threads']['minimap']
    params:
        minimap = config['software']['minimap']
    shell:
        "{params.minimap} -cxasm5 "
            "--paf-no-hit "
            "--cs "
            "-r2k "
            "-y "
            "-t {threads} "
            "{input.ref} {input.fasta} | "
            "gzip -c > {output}"


# Call paf files in var format
# This lets us get unique alignments rather easily for each haplotype
rule call_paf_files:
    input:
        "Diploid_Vars/minimap_contigs{hap}.paf.gz"
    output:
        "Diploid_Vars/minimap_contigs{hap}.var.gz"
    params:
        paftools = config['software']['paftools']
    shell:
        "zcat {input} | "
            "sort -k6,6 -k8,8n | "
            "{params.paftools} call - | "
            "gzip -c > {output}"


# scrape var file for unique alignments andd output as a bed file
rule generate_orthogonal_bed:
    input:
        "Diploid_Vars/minimap_contigs{hap}.var.gz"
    output:
        "Diploid_Vars/minimap_contigs{hap}.bed"
    shell:
        "gzip -dc {input} | "
            "grep ^R | "
            "cut -f2- > {output}"


# run minimap again, this time output a sam file
# these will be used for variant calling
rule run_minimap_sam:
    input:
        fasta = "Assembly/LinkedReads.contig_{hap}.fasta",
        ref = config['params']['ref']
    output:
        "Diploid_Vars/minimap_contigs{hap}.sam.gz"
    params:
        minimap = config['software']['minimap']
    threads:
        config['threads']['minimap']
    shell:
        "{params.minimap} -axasm5 "
            "--paf-no-hit "
            "--cs "
            "-r2k "
            "-t {threads} "
            "{input.ref} {input.fasta} | "
            "gzip -c > {output}"


# sort the sam file, output a bam file
rule flt_sort_to_bam:
    input:
        "Diploid_Vars/minimap_contigs{hap}.sam.gz"
    output:
        "Diploid_Vars/minimap_contigs{hap}.bam"
    params:
        sam_flt = config['software']['sam_flt']
    shell:
        "{params.sam_flt} {input} | "
            "samtools sort -m4G -@4 -o {output}"


# create a bed file of diploid regions from unique alignments
rule bed_diploid_regions:
    input:
        hap0 = "Diploid_Vars/minimap_contigs0.bed",
        hap1 = "Diploid_Vars/minimap_contigs1.bed"
    output:
        "Diploid_Vars/minimap_contig_diploid.bed"
    params:
        bedtools = config['software']['bedtools']
    shell:
        "{params.bedtools} intersect "
            "-a {input.hap0} "
            "-b {input.hap1} > {output}"


# create a haploid bed file of unique alignments
rule bed_haploid_regions:
    input:
        hap0 = "Diploid_Vars/minimap_contigs0.bed",
        hap1 = "Diploid_Vars/minimap_contigs1.bed"
    output:
        "Diploid_Vars/minimap_contig_haploid.bed"
    params:
        bedtools = config['software']['bedtools']
    shell:
        "{params.bedtools} merge -i "
            "<({params.bedtools} sort -i "
                "<({params.bedtools} subtract "
                    "-a {input.hap0} "
                    "-b {input.hap1}; "
                "{params.bedtools} subtract "
                    "-a {input.hap1} "
                    "-b {input.hap0})) > "
            "{output}"


# create a bed file of regions with missing alignments
rule bed_missing_regions:
    input:
        merged = "Diploid_Vars/minimap_contig_diploid.bed",
        haploid = "Diploid_Vars/minimap_contig_haploid.bed"
    output:
        "Diploid_Vars/minimap_contig_missing.bed"
    params:
        bedtools = config['software']['bedtools'],
        ref = config['params']['ref']
    shell:
        "{params.bedtools} complement -i "
            "<({params.bedtools} sort -g "
                "<(cut -f 1-2 {params.ref}.fai) "
                "-i <(cat {input.merged} {input.haploid})) "
            "-g <(cut -f 1-2 {params.ref}.fai) > "
            "{output}"


# function to determine var calling flags
# if -f is present chrom Y is ignored
def is_female(wildcards):
    if config['params']['is_female']:
        return "-f"
    else:
        return ""


# call diploid variants with htsbox pileup
rule diploid_var_calling:
    input:
        hap0 = "Diploid_Vars/minimap_contigs0.bam",
        hap1 = "Diploid_Vars/minimap_contigs1.bam",
    output:
        "Diploid_Vars/minimap_contigs_pair.vcf.gz"
    params:
        htsbox = config['software']['htsbox'],
        ref = config['params']['ref']
    shell:
        "{params.htsbox} pileup -q5 "
            "-evcf "
            "{params.ref} "
            "{input.hap0} "
            "{input.hap1} | "
        "{params.htsbox} bgzip > {output}"


# generate phased VCF with VCF pair
rule phased_vcf:
    input:
        "Diploid_Vars/minimap_contigs_pair.vcf.gz"
    output:
        "Diploid_Vars/minimap_contigs_phased.vcf.gz"
    params:
        vcf_pair = config['software']['vcf_pair'],
        htsbox = config['software']['htsbox'],
        pair_args = is_female
    shell:
        "{params.vcf_pair} {params.pair_args} "
            "{input} | "
        "awk '{{ if ($0 ~ /^#/ || "
            "($4 ~ /^[ACGT]/ && "
             "$5 ~ /^[ACGT]/ )) print }}' | "
        "{params.htsbox} bgzip > {output}"


# index phased VCF
rule phased_vcf_index:
    input:
        "Diploid_Vars/minimap_contigs_phased.vcf.gz"
    output:
        "Diploid_Vars/minimap_contigs_phased.vcf.gz.tbi"
    params:
        htsbox = config['software']['htsbox']
    shell:
        "tabix -p vcf -f {input}"


# returns the wildcard, either snp or indel
def get_type(wildcards):
    type = wildcards.type
    return type


# split VCF into snps and indels, respectively
rule split_vcf:
    input:
        vcf = "Diploid_Vars/minimap_contigs_phased.vcf.gz",
        index = "Diploid_Vars/minimap_contigs_phased.vcf.gz.tbi"
    output:
        "Diploid_Vars/minimap_contig_phased_{type}.vcf.gz"
    params:
        bcftools = config['software']['bcftools'],
        type = get_type
    shell:
        "{params.bcftools} view -O z "
            "--type {params.type}s "
            "{input.vcf} > "
            "{output}"


# index the split VCFs
rule index_split_vcfs:
    input:
        "Diploid_Vars/minimap_contig_phased_{type}.vcf.gz"
    output:
        "Diploid_Vars/minimap_contig_phased_{type}.vcf.gz.tbi"
    shell:
        "tabix -p vcf -f {input}"


# return the appropriate benchmark file for each wildcard type, snp or indel
def benchmark_file(wildcards):
    type = wildcards.type
    if type == "snp":
        return config['benchmark']['benchmark_snp']
    if type == "indel":
        return config['benchmark']['benchmark_indel']


# benchmark snps and indels
# benchmark for diploid, haploid and missing
rule benchmark_regions_and_types:
    input:
        vcf = "Diploid_Vars/minimap_contig_phased_{type}.vcf.gz",
        index = "Diploid_Vars/minimap_contig_phased_{type}.vcf.gz.tbi",
        bed = "Diploid_Vars/minimap_contig_{region}.bed"
    output:
        "Diploid_Vars/minimap_contig_{type}_{region}_vars/summary.txt"
    params:
        truth_vcf = benchmark_file,
        truth_bed = config['benchmark']['bedfile'],
        ref_sdf = config['benchmark']['ref_sdf'],
        rtg = config['software']['rtg']
    shell:
        "rm -r $(dirname {output}); "
        "{params.rtg} vcfeval -b {params.truth_vcf} "
            "--bed-regions={input.bed} "
            "-c {input.vcf} "
            "-e {params.truth_bed} "
            "-t {params.ref_sdf} "
            "-o $(dirname {output})"


# run benchmarking with squashed ploidy
# this looks for the correct nucleotide variant
# and not necessarily the correct ploidy
rule benchmark_regions_and_types_squash_ploidy:
    input:
        vcf = "Diploid_Vars/minimap_contig_phased_{type}.vcf.gz",
        index = "Diploid_Vars/minimap_contig_phased_{type}.vcf.gz.tbi",
        bed = "Diploid_Vars/minimap_contig_{region}.bed"
    output:
        "Diploid_Vars/minimap_squashed_contig_{type}_{region}_vars/summary.txt"
    params:
        truth_vcf = benchmark_file,
        truth_bed = config['benchmark']['bedfile'],
        ref_sdf = config['benchmark']['ref_sdf'],
        rtg = config['software']['rtg']
    shell:
        "rm -r $(dirname {output}); "
        "{params.rtg} vcfeval -b {params.truth_vcf} "
            "--squash-ploidy "
            "--bed-regions={input.bed} "
            "-c {input.vcf} "
            "-e {params.truth_bed} "
            "-t {params.ref_sdf} "
            "-o $(dirname {output})"


# scaffold contigs using python script
rule scaffold_contigs:
    input:
        "Assembly/LinkedReads.contig_{hap}.fasta"
    output:
        "Scaffolds/LinkedReads.scaffolds_{hap}.fa"
    params:
        scaffold = config['software']['scaffold']
    shell:
        "{params.scaffold} {input} > {output}"


# run quast with the scaffolds
rule run_quast_scaffolds:
    input:
        expand("Scaffolds/LinkedReads.scaffolds_{hap}.fa", hap=[0,1])
    output:
        "Scaffolds/Quast_Scaffolds/quast.log"
    threads:
        config['threads']['quast']
    params:
        quast = config['software']['quast'],
        min_contig = config['params']['min_contig'],
        ref = config['params']['ref']
    shell:
        "{params.quast} {input} --eukaryote "
            "-R {params.ref} "
            "--threads {threads} "
            "--no-snps "
            "--fragmented "
            "--min-contig {params.min_contig} "
            "-o $(dirname {output}) "
            "-s "


# generate paf file of scaffolds
rule run_minimap_paf_scaffold:
    input:
        fasta = "Scaffolds/LinkedReads.scaffolds_{hap}.fa",
        ref = config['params']['ref']
    output:
        "Scaffolds/Diploid_Vars/minimap_scaffolds{hap}.paf.gz"
    threads:
        config['threads']['minimap']
    params:
        minimap = config['software']['minimap']
    shell:
        "{params.minimap} -cxasm5 "
            "--paf-no-hit "
            "--cs "
            "-r2k "
            "-y "
            "-t {threads} "
            "{input.ref} {input.fasta} | "
            "gzip -c > {output}"


# call vars in var format from scaffolds
rule call_paf_files_scaffold:
    input:
        "Scaffolds/Diploid_Vars/minimap_scaffolds{hap}.paf.gz"
    output:
        "Scaffolds/Diploid_Vars/minimap_scaffolds{hap}.var.gz"
    params:
        paftools = config['software']['paftools']
    shell:
        "zcat {input} | "
            "sort -k6,6 -k8,8n | "
            "{params.paftools} call - | "
            "gzip -c > {output}"


# scrape the var files for unique alignments
rule generate_orthogonal_bed_scaffold:
    input:
        "Scaffolds/Diploid_Vars/minimap_scaffolds{hap}.var.gz"
    output:
        "Scaffolds/Diploid_Vars/minimap_scaffolds{hap}.bed"
    shell:
        "gzip -dc {input} | "
            "grep ^R | "
            "cut -f2- > {output}"


# run minimap, this time generate sam files
rule run_minimap_sam_scaffold:
    input:
        fasta = "Scaffolds/LinkedReads.scaffolds_{hap}.fa",
        ref = config['params']['ref']
    output:
        "Scaffolds/Diploid_Vars/minimap_scaffolds{hap}.sam.gz"
    params:
        minimap = config['software']['minimap']
    threads:
        config['threads']['minimap']
    shell:
        "{params.minimap} -axasm5 "
            "--paf-no-hit "
            "--cs "
            "-r2k "
            "-y "
            "-t {threads} "
            "{input.ref} {input.fasta} | "
            "gzip -c > {output}"


# sort sam file and output bam
rule flt_sort_to_bam_scaffold:
    input:
        "Scaffolds/Diploid_Vars/minimap_scaffolds{hap}.sam.gz"
    output:
        "Scaffolds/Diploid_Vars/minimap_scaffolds{hap}.bam"
    params:
        sam_flt = config['software']['sam_flt']
    shell:
        "{params.sam_flt} {input} | "
            "samtools sort -m4G -@4 -o {output}"


# generate diploid bed file of scaffold unique alignments
rule bed_diploid_regions_scaffold:
    input:
        hap0 = "Scaffolds/Diploid_Vars/minimap_scaffolds0.bed",
        hap1 = "Scaffolds/Diploid_Vars/minimap_scaffolds1.bed"
    output:
        "Scaffolds/Diploid_Vars/minimap_scaffold_diploid.bed"
    params:
        bedtools = config['software']['bedtools']
    shell:
        "{params.bedtools} intersect "
            "-a {input.hap0} "
            "-b {input.hap1} > {output}"


# generate haploid bed file of scaffold alignments
rule bed_haploid_regions_scaffold:
    input:
        hap0 = "Scaffolds/Diploid_Vars/minimap_scaffolds0.bed",
        hap1 = "Scaffolds/Diploid_Vars/minimap_scaffolds1.bed"
    output:
        "Scaffolds/Diploid_Vars/minimap_scaffold_haploid.bed"
    params:
        bedtools = config['software']['bedtools']
    shell:
        "{params.bedtools} merge -i "
            "<({params.bedtools} sort -i "
                "<({params.bedtools} subtract "
                    "-a {input.hap0} "
                    "-b {input.hap1}; "
                "{params.bedtools} subtract "
                    "-a {input.hap1} "
                    "-b {input.hap0})) > "
            "{output}"


# generate missing bed file of scaffold alignments
rule bed_missing_regions_scaffold:
    input:
        merged = "Scaffolds/Diploid_Vars/minimap_scaffold_diploid.bed",
        haploid = "Scaffolds/Diploid_Vars/minimap_scaffold_haploid.bed"
    output:
        "Scaffolds/Diploid_Vars/minimap_scaffold_missing.bed"
    params:
        bedtools = config['software']['bedtools'],
        ref = config['params']['ref']
    shell:
        "{params.bedtools} complement -i "
            "<({params.bedtools} sort -g "
                "<(cut -f 1-2 {params.ref}.fai) "
                "-i <(cat {input.merged} {input.haploid})) "
            "-g <(cut -f 1-2 {params.ref}.fai) > "
            "{output}"


# call vars using htsbox pileup
rule diploid_var_calling_scaffold:
    input:
        hap0 = "Scaffolds/Diploid_Vars/minimap_scaffolds0.bam",
        hap1 = "Scaffolds/Diploid_Vars/minimap_scaffolds1.bam"
    output:
        "Scaffolds/Diploid_Vars/minimap_scaffolds_pair.vcf.gz"
    params:
        htsbox = config['software']['htsbox'],
        ref = config['params']['ref']
    shell:
        "{params.htsbox} pileup -q5 "
            "-evcf "
            "{params.ref} "
            "{input.hap0} "
            "{input.hap1} | "
        "{params.htsbox} bgzip > {output}"


# phase VCF with vcf_pair
rule phased_vcf_scaffold:
    input:
        "Scaffolds/Diploid_Vars/minimap_scaffolds_pair.vcf.gz"
    output:
        "Scaffolds/Diploid_Vars/minimap_scaffolds_phased.vcf.gz"
    params:
        vcf_pair = config['software']['vcf_pair'],
        htsbox = config['software']['htsbox'],
        pair_args = is_female
    shell:
        "{params.vcf_pair} {params.pair_args} "
            "{input} | "
        "awk '{{ if ($0 ~ /^#/ || "
            "($4 ~ /^[ACGT]/ && "
             "$5 ~ /^[ACGT]/ )) print }}' | "
        "{params.htsbox} bgzip > {output}"


# index the phased VCF
rule phased_vcf_index_scaffold:
    input:
        "Scaffolds/Diploid_Vars/minimap_scaffolds_phased.vcf.gz"
    output:
        "Scaffolds/Diploid_Vars/minimap_scaffolds_phased.vcf.gz.tbi"
    params:
        htsbox = config['software']['htsbox']
    shell:
        "tabix -p vcf -f {input}"


# split scaffold VCF into SNPs and InDels
rule split_vcf_scaffold:
    input:
        vcf = "Scaffolds/Diploid_Vars/minimap_scaffolds_phased.vcf.gz",
        index = "Scaffolds/Diploid_Vars/minimap_scaffolds_phased.vcf.gz.tbi"
    output:
        "Scaffolds/Diploid_Vars/minimap_scaffold_phased_{type}.vcf.gz"
    params:
        bcftools = config['software']['bcftools'],
        type = get_type
    shell:
        "{params.bcftools} view -O z "
            "--type {params.type}s "
            "{input.vcf} > "
            "{output}"


# index the split VCFs
rule index_split_vcfs_scaffold:
    input:
        "Scaffolds/Diploid_Vars/minimap_scaffold_phased_{type}.vcf.gz"
    output:
        "Scaffolds/Diploid_Vars/minimap_scaffold_phased_{type}.vcf.gz.tbi"
    shell:
        "tabix -p vcf -f {input}"


# run benchmarking for SNPs and InDels
# for haploid, diploid and missing regions
rule benchmark_regions_and_types_scaffold:
    input:
        vcf = "Scaffolds/Diploid_Vars/minimap_scaffold_phased_{type}.vcf.gz",
        index = "Scaffolds/Diploid_Vars/minimap_scaffold_phased_{type}.vcf.gz.tbi",
        bed = "Scaffolds/Diploid_Vars/minimap_scaffold_{region}.bed"
    output:
        "Scaffolds/Diploid_Vars/minimap_scaffold_{type}_{region}_vars/summary.txt"
    params:
        truth_vcf = benchmark_file,
        truth_bed = config['benchmark']['bedfile'],
        ref_sdf = config['benchmark']['ref_sdf'],
        rtg = config['software']['rtg']
    shell:
        "rm -r $(dirname {output}); "
        "{params.rtg} vcfeval -b {params.truth_vcf} "
            "--bed-regions={input.bed} "
            "-c {input.vcf} "
            "-e {params.truth_bed} "
            "-t {params.ref_sdf} "
            "-o $(dirname {output})"


# run again with ploidy squashed
rule benchmark_regions_and_types_squash_ploidy_scaffold:
    input:
        vcf = "Scaffolds/Diploid_Vars/minimap_scaffold_phased_{type}.vcf.gz",
        index = "Scaffolds/Diploid_Vars/minimap_scaffold_phased_{type}.vcf.gz.tbi",
        bed = "Scaffolds/Diploid_Vars/minimap_scaffold_{region}.bed"
    output:
        "Scaffolds/Diploid_Vars/minimap_squashed_scaffold_{type}_{region}_vars/summary.txt"
    params:
        truth_vcf = benchmark_file,
        truth_bed = config['benchmark']['bedfile'],
        ref_sdf = config['benchmark']['ref_sdf'],
        rtg = config['software']['rtg']
    shell:
        "rm -r $(dirname {output}); "
        "{params.rtg} vcfeval -b {params.truth_vcf} "
            "--squash-ploidy "
            "--bed-regions={input.bed} "
            "-c {input.vcf} "
            "-e {params.truth_bed} "
            "-t {params.ref_sdf} "
            "-o $(dirname {output})"
