import os

rule fastqc:
    input:
        r1 = "resources/{sample}_R1.fastq.gz",
        r2 = "resources/{sample}_R2.fastq.gz"
    output:
        "results/fastqc/{sample}_R1_fastqc.zip",
        "results/fastqc/{sample}_R2_fastqc.zip",
        "results/fastqc/{sample}_R1_fastqc.html",
        "results/fastqc/{sample}_R2_fastqc.html"
    log:
        "results/logs/{sample}.fastqc.log"
    conda:
        "../envs/merip.yaml"
    priority: 1
    params:
        out_dir=lambda wildcards, output: os.path.split(output[0])[0]
    benchmark:
        "results/benchmarks/{sample}.fastqc.benchmark.txt"
    threads: 2
    shell:
        """
        fastqc -t {threads} {input.r1} {input.r2} \
            -o {params.out_dir} > {log} 2>&1
        """