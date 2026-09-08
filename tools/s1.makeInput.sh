#!/usr/bin/env bash
# yhbatch -N 1 -n 64 -p deimos
# yhbatch -p e9 -n 1
pheno=ND_13_ALL.txt
outprefix=ND_13_ALL
num=80

cat $pheno |sed '1d'|cut -f1|head -n $num |while read id;do echo "/XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/20260807_RRBS/hsa/output/01.fq2methRatio/$id/03.methylation/$id.CG_5x.bedgraph.gz";done > g1.list

cat $pheno |sed '1d'|cut -f1|sed "1,${80}d"|while read id;do echo "/XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/20260807_RRBS/hsa/output/01.fq2methRatio/$id/03.methylation/$id.CG_5x.bedgraph.gz";done > g2.list

# metilene_input=/XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_1/BIGDATA2/gzfezx_shhli_1/softwares/metilene_v0.2-8/metilene_input_plus.pl
mkinputScript=/XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/script/rrbs/tools/metilene_input_plus.py
bedtools=/XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_1/BIGDATA2/gzfezx_shhli_1/softwares/bedtools2/bin/bedtools

time python $mkinputScript --in1_list g1.list --in2_list g2.list --out ${outprefix}_input_R01.txt -b $bedtools -t 13 --ratio 0.1 && rm g1.list g2.list

# echo -e "chrom pos "$(less $pheno |sed '1d'|cut -d"," -f1|tr '\n' ' ')> header.txt
# sed -i 's/ /\t/g' header.txt 
# sed '1d' ${outprefix}_input_R01.txt  | cat header.txt -  > ${outprefix}_input_R01.bigcsID.txt && rm header.txt g1.list g2.list
