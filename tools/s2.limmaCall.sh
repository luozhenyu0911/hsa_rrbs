#!/usr/bin/env bash
# yhbatch -N 1 -n 64 -p deimos
# yhbatch -p e9 -n 1

# 生成xx个case和xx个control, 共336个样本(必须是按照case和control的顺序排列的)
str="$(printf 'case,%.0s' {1..80})$(printf 'control,%.0s' {1..179})"
str="${str%,}"

inputfile=ND_13_ALL_input_R01.txt
outputDir=ND_13_ALL
outputfilePrefix=case_vs_control_07
min_coverage=0.7

limma_diff=/XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/script/rrbs/tools/limma_diff.R
cellprops=/XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/20260807_RRBS/hsa/output/02.s856CBcorrect/s856CB_methRatio.CellRatio.CB.uniLIFE7.txt

time Rscript $limma_diff \
  -i $inputfile \
  -g $str \
  -o $outputDir \
  --prefix $outputfilePrefix \
  --dmp_fdr 0.05 \
  --dmp_delta 0.1 \
  --dmr_fdr 0.05 \
  --dmr_meandiff 0.1 \
  --lambda 1000 \
  --C 2 \
  --coverage_by group \
  --min_coverage $min_coverage \
  --impute group \
  --cell_props $cellprops && \
  echo "done"
