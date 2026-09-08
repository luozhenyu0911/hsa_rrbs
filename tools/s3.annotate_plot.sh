#!/bin/bash
# yhbatch -N 1 -n 64 -p deimos
# yhbatch -p e9 -n 1
echo "start time: $(date '+%H:%M:%S')"
start_time=$(date +%s)  # record start time

source /HOME/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/miniconda3/bin/activate manta2 && echo "has activated manta2 environment of gzfezx_shhli_2" && \

# setting config file 
CONFIG_FILE="./config.conf"
# source config file
source "$CONFIG_FILE"

files=(
    "$calling_dir/$name.DMC_pval.0.05.out"
    "$calling_dir/$name.DMR_qval.0.05.out"
)

for file in "${files[@]}"; do
    if [ ! -s "$file" ]; then
        echo "the file is missing or empty: $file"
        exit 1
    fi
done


/HOME/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/miniconda3/bin/python \
                /XYFS02/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_1/Projects/20260514_S256rrbs/bin/calbedmeth_v5.py \
                -m $calling_dir/$name.inputf.txt \
                -b $calling_dir/$name.DMR_qval.0.05.out \
                -o $calling_dir/$name.DMR_qval.0.05.sampleMet.bed && \

Rscript ../bin/plot1.annotate_DMR.R $name > plot1.annotate_DMR.log
Rscript ../bin/plot2.gsea_DMR.R $name > plot2.gsea_DMR.log 
Rscript ../bin/plot3.heatmap_DMR.R $name $pheno_file > plot3.heatmap_DMR.log
Rscript ../bin/plot4.manhattan_DMC.R $name > plot4.manhattan_DMC.log

mv *plot*pdf 02.anno_plot/Out_Manhtn

end_time=$(date +%s)
echo "end time: $(date '+%H:%M:%S')"
duration=$((end_time - start_time))
hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))
seconds=$((duration % 60))
printf "all processes time: %02d:%02d:%02d\n" $hours $minutes $seconds

echo "##############################   finished  ############################################"












