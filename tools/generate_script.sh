#!/usr/bin/env bash
# yhbatch -N 1 -n 64 -p deimos
# yhbatch -p e9 -n 1
python generate_script.py  \
-i input_sh/NIPT_s1k.list -c config.yaml \
-r run_snakemake.sh -p /XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/20260506_NIPT_viromes
