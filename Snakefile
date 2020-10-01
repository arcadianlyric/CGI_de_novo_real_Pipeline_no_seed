# Snakemake Pipeline for real de novo datasets without a seed library
# Start from absolute scratch
# So need to make the assembly using the newest software (config file)
# prepend commands with export Sentieon stuff.

configfile: "config.yaml"
# These export the sentieon environment variables for all commands
shell.prefix('export SENTIEON_INSTALL={}; export SENTIEON_LICENSE={}; export SENTIEON_TMPDIR={};'.format(config['params']['sentieon_install'], config['params']['sentieon_license'], config['params']['tmp_dir']))

# Rule for run all, could be shortened with a function
# These are the targets we want to make
rule run_all:
    input:
        "Quast_Contigs/quast.log",
        "hap0_phase_n50.txt",
        "hap1_phase_n50.txt",
        "Assembly/LinkedReads.seed_0.fasta",
        "Assembly/LinkedReads.seed_1.fasta",
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


# rule for read trimming using atropos
rule atropos:
    input:
        expand("{sample}", sample=config['samples']['fastq'])
    output:
        "atropos_interleaved.fq.gz"
    params:
        adapt = config['atropos']['adapter'],
        adapt2 = config['atropos']['adapter2'],
        front = config['atropos']['front'],
        front2 = config['atropos']['front2'],
        min_len = config['atropos']['min_len'],
        vir_env = config['atropos']['vir_env']
    threads:
        config['threads']['atropos']
    run:
        virtual_env = str(params.vir_env) + ';'
        command = ['source', virtual_env, 'atropos', 'trim', '-T', str(threads)
                   '-m', str(params.min_len), '-pe1', input[0], '-pe2', input[1],
                   '-L', '/dev/stdout']

        # if adapters are included, add them to the command
        if params.adapt:
            command.extend(['--adapter', params.adapt])

        if params.adapt2:
            command.extend(['--adapter2', params.adapt2])

        if params.front:
            command.extend(['--front', params.front])

        if params.front2:
            command.extend(['--front2', params.front2])

        command.extend(['|', 'gzip', '-c', '>', output[0]])
        # run the command
        shell(" ".join(command))


# create unitigs
rule run_bcalm:
    input:
        "atropos_interleaved.fq.gz"
    output:
        "Assembly/reads.unitigs.fa"
    params:
        kmer_size = config['params']['ksize'],
        a_min = config['params']['amin']
    threads:
        config['threads']['assembly']
    log:
        "Assembly/bcalm.err"
    shell:
        "$SENTIEON_INSTALL/bin/sentieon bcalm "
            "-kmer-size {params.kmer_size} "
            "-abundance-min {params.a_min} "
            "-in {input} "
            "-out Assembly/reads "
            "-nb-cores {threads} > {log} 2>&1"


# index unitigs with samtools
rule index_unitigs_samtools:
    input:
        "Assembly/reads.unitigs.fa"
    output:
        "Assembly/reads.unitigs.fa.fai"
    shell:
        "samtools faidx {input}"


# index unitigs with bwa
rule index_unitigs_bwa:
    input:
        "Assembly/reads.unitigs.fa"
    output:
        "Assembly/reads.unitigs.fa.bwt"
    log:
        "Assembly/bwa.index.err"
    shell:
        "$SENTIEON_INSTALL/bin/sentieon bwa index {input} > {log} 2>&1"


# assign reads to unitigs
# -C will append fq comments to the outputs, -p handles interleaved reads
rule assign_reads_to_unitigs:
    input:
        fasta = "Assembly/reads.unitigs.fa",
        fai = "Assembly/reads.unitigs.fa.fai",
        bwt = "Assembly/reads.unitigs.fa.bwt",
        reads = "atropos_interleaved.fq.gz"
    output:
        "Assembly/aligned.bam"
    params:
        chunk = config['params']['chunk'],
        ksize = config['params']['ksize'],
        readgroup = r'@RG\\tID:{0}\\tSM:{0}\\tPL:{1}'.format(config['samples']['sample'],
                                                          config['params']['platform'])
    log:
        bwa = "Assembly/bwa.reads.err",
        group = "Assembly/group.reads.err"
    threads:
        config['threads']['assembly']
    shell:
        "$SENTIEON_INSTALL/bin/sentieon bwa mem "
            "-R {params.readgroup} "
            "-t {threads} "
            "-K {params.chunk} "
            "-k {params.ksize} "
            "-x unitig "
            "-p {input.fasta} "
            "{input.reads} 2> {log.bwa} | "
        "extract | "
        "ASAN_OPTIONS=detect_leaks=0 $SENTIEON_INSTALL/bin/sentieon driver "
            "-t {threads} "
            "--algo LinkedReadsDeNovo group "
            "{output} - > {log.group} 2>&1"


def get_seed_sizes():
    seed_sizes = [str(x) for x in config['params']['seed_sizes']]
    return ",".join(seed_sizes)


# perform de novo assembly
rule de_novo_assembly:
    input:
        bam = "Assembly/aligned.bam",
        fasta = "Assembly/reads.unitigs.fa"
    output:
        contigs = expand("Assembly/LinkedReads.seed_{hap}.fasta", hap=[0,1]),
        phase_contigs = expand("Assembly/LinkedReads.seed.phase_{hap}.fasta", hap=[0,1])
    params:
        ksize = config['params']['ksize'],
        trace_size = config['params']['tracegraph_size'],
        read_len = config['params']['read_len'],
        seed_sizes = get_seed_sizes()
    threads:
        config['threads']['assembly']
    log:
        "Assembly/LinkedReadsSeed.log"
    shell:
        "$SENTIEON_INSTALL/bin/sentieon driver -t {threads} "
            "--algo LinkedReadsDeNovo trace "
            "-i {input.bam} "
            "--bcalm_untig_graph {input.fasta} "
            "--untig_kmer_size {params.ksize} "
            "--contig_output Assembly/LinkedRead.seed "
            "--seed_sizes {params.seed_sizes} "
            "--max_tracegraph_size {params.trace_size} "
            "--read_len {params.read_len} > {log} 2>&1"


# run quast on contig assemblies
rule run_quast_contigs:
    input:
        expand("Assembly/LinkedReads.seed_{hap}.fasta", hap=[0,1])
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
        "Assembly/LinkedReads.seed.phase_{hap}.fasta"
    output:
        "hap{hap}_phase_n50.txt"
    params:
        phase_script = config['software']['phase_script']
    shell:
        "{params.phase_script} {input} > {output}"


# Generate paf files using mimimap
rule run_minimap_paf:
    input:
        fasta = "Assembly/LinkedReads.seed_{hap}.fasta",
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
        fasta = "Assembly/LinkedReads.seed_{hap}.fasta",
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
        "Assembly/LinkedReads.seed_{hap}.fasta"
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
