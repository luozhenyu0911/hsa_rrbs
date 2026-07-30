
rule summary_conversion_rate:
    input:
        "04.metrics/{id}_lambda_PE_report.txt"
    output:
        "04.metrics/{id}_conversion_rate.txt"
    shell:
        """
        time python3 /XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/script/rrbs/tools/calcu_conversion_rate.py {input} {output}
        """