# 6. macs2 call peak
rule macs2:
    input:
        Ip="results/hisat2/{sample}-IP.bam",
        Input="results/hisat2/{sample}-INPUT.bam"
    output:
        "results/macs2/{sample}/m6a_peaks.narrowPeak"
    params:
        outdir= "results/macs2/{sample}"
    log:
        "results/logs/{sample}.macs2.log"
    conda:
        "../envs/macs2.yaml"
    benchmark:
        "results/benchmarks/{sample}.macs2.benchmark.txt"
    shell:
        """
        macs2 callpeak \
            -t {input.Ip} \
            -c {input.Input} \
            --nomodel -f BAMPE \
            -g mm -n m6a -q 0.05 -B\
            --outdir {params.outdir}> {log} 2>&1
        """