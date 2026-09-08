#!/bin/bash
# yhbatch -N 1 -n 64 -p deimos
# # yhbatch -p e9 -n 1
source /XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/miniconda3/bin/activate rrbs && echo "rrbs environment activated." && \

pwd=/XYFS02/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/Project/test_snakemake
smk=/XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/script/rrbs/smkfiles/run.all.smk
threads=8

time ~/BIGDATA2/gzfezx_shhli_2/miniconda3/envs/cnvnator/bin/snakemake -d $pwd -c $threads -pk --configfile ${pwd}/config.yaml --default-resources tmpdir="\"${pwd}\"" -s ${smk} --rerun-triggers mtime > ${pwd}/snakemake.err.txt 2>&1 

