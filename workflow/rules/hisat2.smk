# 4. hisat2比对
rule hisat2:
    input:
        # r1 = "results/trim/{sample}.trim_R1.fastq.gz",
        # r2 = "results/trim/{sample}.trim_R2.fastq.gz",
        r1 = "results/cleanData/{sample}_cleanData.fastq.1.gz",
        r2 = "results/cleanData/{sample}_cleanData.fastq.2.gz"
    output:
        temp("results/hisat2/{sample}.sam")
    log:
        "results/logs/{sample}.hisat2.log"
    conda:
        "../envs/merip.yaml"
    threads: workflow.cores * 0.2
    benchmark:
        "results/benchmarks/{sample}.hisat2.benchmark.txt"
    params:
        index = config["index"]
    shell:
        """
        hisat2 -p {threads} -N 1 --dta -x {params} \
            -1 {input.r1} \
            -2 {input.r2} \
            -S  {output}> {log} 2>&1
        """
rule sam2bam:
    input:
        "results/hisat2/{sample}.sam"
    output:
        protected("results/hisat2/{sample}.bam")
    log:
        "results/logs/{sample}.sam2bam.log"
    conda:
        "../envs/merip.yaml"
    threads: workflow.cores * 0.2
    benchmark:
        "results/benchmarks/{sample}.sam2bam.benchmark.txt"
    shell:
        """
        samtools view {input} -q {threads} -h -S -b -o {output}> {log} 2>&1
        """
ruleorder: sorted_bam > index_bam > bamCoverage
rule sorted_bam:
    input:
        "results/hisat2/{sample}.bam"
    output:
        "results/hisat2/{sample}.sorted.bam"
    log:
        "results/logs/{sample}.sorted_bam.log"
    conda:
        "../envs/merip.yaml"
    threads: workflow.cores * 0.2
    benchmark:
        "results/benchmarks/{sample}.sorted_bam.benchmark.txt"
    shell:
        """
        samtools sort {input} -o {output} -@ {threads} > {log} 2>&1
        """
rule index_bam:
    input:
        "results/hisat2/{sample}.sorted.bam"
    output:
        "results/hisat2/{sample}.sorted.bam.bai"
    log:
        "results/logs/{sample}.index_bam.log"
    conda:
        "../envs/merip.yaml"
    threads: workflow.cores * 0.2
    priority: 4
    benchmark:
        "results/benchmarks/{sample}.index_bam.benchmark.txt"
    shell:
        """
        samtools index {input} -@ {threads} > {log} 2>&1
        """
rule bamCoverage:
    input:
        "results/hisat2/{sample}.sorted.bam"
    output:
        "results/bigwig/{sample}.bigwig"
    log:
        "results/logs/{sample}.bamCoverage.log"
    conda:
        "../envs/deeptools.yaml"
    threads: workflow.cores * 0.2
    benchmark:
        "results/benchmarks/{sample}.bamCoverage.benchmark.txt"
    shell:
        """
        bamCoverage -p {threads} --bam {input} --normalizeUsing RPKM --binSize 10 --outFileName {output} > {log} 2>&1
        """