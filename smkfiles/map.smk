
from pathlib import Path

rule link_clean_fq1:
    input:
        config['samples']['raw_fq1']
    output:
        "data/{id}.R1.fq.gz"
    run:
        shell("ln -s {input} {output}")

rule link_clean_fq2:
    input:
        config['samples']['raw_fq2']
    output:
        "data/{id}.R2.fq.gz"
    run:
        shell("ln -s {input} {output}")

rule qc_trim:
    input:
        fq1 = "data/{id}.R1.fq.gz",
        fq2 = "data/{id}.R2.fq.gz"
    output:
        fq1 = "01.qc/{id}_val_1.fq.gz",
        fq2 = "01.qc/{id}_val_2.fq.gz",
        R1_html = "01.qc/{id}_val_1_fastqc.html",
        R2_html = "01.qc/{id}_val_2_fastqc.html"
    threads:
        config['threads']
    params:
        dir = "01.qc/"
    shell:
        """
        time trim_galore \
            --paired \
            --rrbs \
            --quality 20 \
            --stringency 3 \
            --length 20 \
            --cores {threads} \
            --fastqc_args '--threads {threads}' \
            --basename {wildcards.id} \
            -o {params.dir} \
            {input.fq1} {input.fq2}
        """

rule qc_fastqc:
    input:
        fq1 = "data/{id}.R1.fq.gz",
        fq2 = "data/{id}.R2.fq.gz"
    output:
        R1_html = "01.qc/{id}.R1_fastqc.html",
        R2_html = "01.qc/{id}.R2_fastqc.html"
    threads:
        config['threads']
    shell:
        "time fastqc -t {threads} -o 01.qc/ {input.fq1} {input.fq2}"

rule qc_multiqc:
    input:
        R1_html = "01.qc/{id}_val_1_fastqc.html",
        R2_html = "01.qc/{id}_val_2_fastqc.html"
    output:
        multiqc_html = "01.qc/{id}.multiqc.html"
    threads:
        config['threads']
    params:
        filename = '{id}.multiqc'
    shell:
        """
        time multiqc 01.qc/ -o 01.qc/ -n {params.filename} && rm -f 01.qc/{wildcards.id}_val_*.fq.gz && rm -rf 01.qc/{wildcards.id}.multiqc_data && rm -f 01.qc/{wildcards.id}*_fastqc.zip 
        """

rule mapping2lambda:
    input:
        # fq1 = "01.qc/{id}_val_1.fq.gz",
        # fq2 = "01.qc/{id}_val_2.fq.gz",
        fq1 = "data/{id}.R1.fq.gz",
        fq2 = "data/{id}.R2.fq.gz",
        lambda_REF = config['params']['ref_lambda']  # 确保此路径存在且是索引文件夹
    output:
        report = "04.metrics/{id}_lambda_PE_report.txt"
        # bam = "04.metrics/{id}_lambda_pe.bam"
    threads:
        config['threads']
    params:
        output_dir = "04.metrics/",
        prefix = "{id}_lambda"
    shell:
        """
        time seqtk sample -s100 {input.fq1} 5000000 |gzip > 04.metrics/{wildcards.id}.R1.fq.gz & \
        time seqtk sample -s100 {input.fq2} 5000000 |gzip > 04.metrics/{wildcards.id}.R2.fq.gz & \
        wait && \
        time bismark --bowtie2 -p {threads} --output_dir {params.output_dir} \
            --basename {params.prefix} --temp_dir {params.output_dir} \
            {input.lambda_REF} \
            -1 04.metrics/{wildcards.id}.R1.fq.gz -2 04.metrics/{wildcards.id}.R2.fq.gz \
        && rm -f 04.metrics/{wildcards.id}_lambda_pe.bam && rm -f 04.metrics/{wildcards.id}.R*.fq.gz
        """

rule mapping2hg38:
    input:
        fq1 = "01.qc/{id}_val_1.fq.gz",
        fq2 = "01.qc/{id}_val_2.fq.gz",
        REF = config['params']['ref_hg38']
    output:
        report = "02.bismark_bt2/{id}_PE_report.txt",
        bam = "02.bismark_bt2/{id}_pe.bam"
    threads:
        config['threads']
    params:
        output_dir = "02.bismark_bt2/"
    shell:
        """
        time bismark --bowtie2 -p {threads} --output_dir {params.output_dir} \
            --basename {wildcards.id} --temp_dir {params.output_dir} \
            {input.REF} \
            -1 {input.fq1} -2 {input.fq2}
        """

rule extract_methylation:
    input:
        bam = "02.bismark_bt2/{id}_pe.bam",
        REF = config['params']['ref_hg38']
    output:
        # report = "03.methylation/{id}_pe_splitting_report.txt",
        bedgraph = "03.methylation/{id}_pe.bedGraph.gz",
        cov = "03.methylation/{id}_pe.bismark.cov.gz"
    threads:
        config['threads']
    shell:
        """
        time bismark_methylation_extractor --gzip --bedGraph --no_overlap --comprehensive \
        --parallel {threads} \
        --genome_folder {input.REF} \
        -o 03.methylation \
        {input.bam}
        # && rm -f {input.bam}
        """

rule extract_5x_bedGraph:
    input:
        cov = "03.methylation/{id}_pe.bismark.cov.gz"
    output:
        bedgraph = "03.methylation/{id}.CG_5x.bedgraph.gz",
        stat_file = "04.metrics/{id}.CG_depth.stat.txt"
    params:
        tools = config['params']['tools']
    shell:
        """
        time python {params.tools}/cov2bedGraph.py -r {input.cov} -b {output.bedgraph} -s {output.stat_file} && rm -f 03.methylation/*_pe.txt.gz
        """

