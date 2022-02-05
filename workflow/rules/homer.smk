
# # samtools faidx genome.fa得到sizes文件
# shuffleBed -i $OUTPUT/ok.bed -g /home/shpc_93816f5f84/projects/TRF-m6A/data/GRCm38.p6.chrom.sizes  > $OUTPUT/shuffle.bed

# findMotifsGenome.pl $OUTPUT/ok.bed \
#     /home/shpc_93816f5f84/ref/refdata-gex-mm10-2020-A/fasta/genome.fa \
#     $OUTPUT/TOTAL -size given -p 20 -len 5,6,7,8 -rna -chopify \
#     -norevopp -cache 1000 -bg $OUTPUT/shuffle.bed

# 7. homer 找motif
rule get_ok_bed:
    input:
        "results/macs2/{sample}/m6a_peaks.narrowPeak"
    output:
        "results/macs2/{sample}/ok.bed"
    shell:
        """
        cat {input} |awk -v OFS="\t" '{{print $1,$2,$3,"+"}}'> {output}
        """
rule get_shuffle_bed:
    input:
        "results/macs2/{sample}/ok.bed"
    output:
        "results/macs2/{sample}/shuffle.bed"
    conda:
        "../envs/bedtools.yaml"
    params:
        chrom_sizes= config["chrom_sizes"]
    shell:
        """
        shuffleBed -i {input} -g {params.chrom_sizes}  > {output}
        """

rule findMotifsGenome:
    input:
        shuffle= "results/macs2/{sample}/shuffle.bed",
        ok = "results/macs2/{sample}/ok.bed"
    output:
        "results/macs2/{sample}/TOTAL/homerResults.html"
    params:
        outdir= "results/macs2/{sample}/TOTAL",
        genome = config["genome"]
    log:
        "results/logs/{sample}.findMotifsGenome.log"
    conda:
        "../envs/homer.yaml"
    benchmark:
        "results/benchmarks/{sample}.findMotifsGenome.benchmark.txt"
    threads: workflow.cores
    shell:
        """
        findMotifsGenome.pl {input.ok} \
            {params.genome} {params.outdir} \
            -size 200 -p {threads} -len 6,8,10 \
             -bg {input.shuffle}  > {log} 2>&1
        """