
rule summary_conversion_rate:
    input:
        "04.metrics/{id}_lambda_PE_report.txt"
    output:
        "04.metrics/{id}_conversion_rate.txt"
    shell:
        """
        time python3 /XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/script/rrbs/tools/calcu_conversion_rate.py {input} {output}
        """

rule bam2cram:
    input:
        bam = "02.bismark_bt2/{id}_pe.bam",
        REF = config['params']['ref_hg38']
    output:
       "02.bismark_bt2/{id}_pe.cram"
    threads:
        config['threads']
    
    shell:
        """
        samtools view -@ {threads} -C -T /XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/database/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa -o {output} {input.bam} && rm -f {input.bam}
        """