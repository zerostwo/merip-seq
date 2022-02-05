# 移除rRNA
rule rm_rRNA:
    input:
        r1 = "results/trim/{sample}.trim_R1.fastq.gz",
        r2 = "results/trim/{sample}.trim_R2.fastq.gz"
    output:
        "results/cleanData/{sample}_cleanData.fastq.1.gz",
        "results/cleanData/{sample}_cleanData.fastq.2.gz"
    log:
        "results/logs/{sample}.rm_rRNA.log"
    threads: workflow.cores * 0.2
    benchmark:
        "results/benchmarks/{sample}.rm_rRNA.benchmark.txt"
    params:
        rRNA_index = config["rRNA_index"],
        output_file = "results/cleanData/{sample}_cleanData.fastq.gz",
        sam = "results/cleanData/{sample}.sam"
    shell:
        """
        ~/soft/bowtie2-2.4.4-linux-x86_64/bowtie2 -x {params.rRNA_index} \
            --un-conc-gz {params.output_file} \
            -1 {input.r1} \
            -2 {input.r2} \
            -p {threads} -S {params.sam}; rm {params.sam} > {log} 2>&1
        """