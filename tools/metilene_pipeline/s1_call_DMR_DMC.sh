#!/bin/bash
# yhbatch -N 1 -n 64 -p deimos
# yhbatch -p e9 -n 1
echo "start time: $(date '+%H:%M:%S')"
start_time=$(date +%s)  # record start time

source /HOME/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/miniconda3/bin/activate manta2 && echo "has activated manta2 environment of gzfezx_shhli_2"

# setting config file 
CONFIG_FILE="./config.conf"
# source config file
source "$CONFIG_FILE"

# check input files
for file in "$CONFIG_FILE" "$pheno_file" "$inputfile"; do
    if [ ! -f "$file" ]; then
        echo "Error: Required file '$file' not found"
        exit 1
    fi
done
echo -e "  $CONFIG_FILE\n  $pheno_file\n  $inputfile\nAll three of the above files exist."

echo ""
echo "##############################   preparing input files  #####################################"
echo "###  step1:according to your pheno.csv,generate the grouped matrix file"
dos2unix $pheno_file

# check whether the first value of the column in pheno_file matches the name
if [ "$(awk -F"," -v col="$column" 'NR==1{print $col; exit}' $pheno_file)" = $name ]; then
    echo "Match: First row column $column equals '$name' in $pheno_file"
else
    echo "Error: First row column $column does not equal '$name' in $pheno_file, please check your setting in config.conf"
    exit 1
fi

# Check whether the sample name order in pheno_file matches that of the original pheno file.
#!/bin/bash
if ! diff -q <(awk -F"," '{print $2}' "$pheno_file") <(awk -F"," '{print $2}' "$original_pheno_file") >/dev/null; then
    echo "Error: $pheno_file File order does not match the original file $original_pheno_file"
    echo "Please check: $pheno_file and $original_pheno_file"
    exit 1
    else
        echo "Match: $pheno_file File order matches the original file $original_pheno_file"
        echo ""
fi

mkdir -p $calling_dir
cat $pheno_file |sed '1d'|awk -v col="$column" -F',' '{print $col}'|sed '1i chrom\npos'| tr '\n' '\t'| sed 's/\t$/\n/' > $calling_dir/$name.header.txt

sed '1d' $inputfile | cat $calling_dir/$name.header.txt -  > $calling_dir/$name.inputf.txt

echo "##############################   for DMR/DMC analysis  ############################################"
echo "###  step2:use metilene to test the difference between two groups, results are generate in outfile directory"

# metilene -a g1 -b g2 -t $threads $calling_dir/$name.inputf.txt | sort -k1,1V -k2,2n  > $calling_dir/$name.DMR.bed && \ 
# metilene -a g1 -b g2 -f 3 -t $threads $calling_dir/$name.inputf.txt | sort -k1,1V -k2,2n  > $calling_dir/$name.DMC.bed && \

# 检查 sort 版本
# sort --version | head -1
# GNU sort 8.x 支持 --parallel

metilene -a g1 -b g2 -t $threads "$calling_dir/$name.inputf.txt" > "$calling_dir/$name.DMR.raw.txt" 2>"$calling_dir/$name.DMR.log" && \
sort --temporary-directory="./" --parallel=$threads -k1,1V -k2,2n "$calling_dir/$name.DMR.raw.txt" > "$calling_dir/$name.DMR.bed" && \
perl $tools/metilene_output.pl -q $calling_dir/$name.DMR.bed \
                               -o $calling_dir/$name.DMR \
                               -p 0.05 -c 5 -d 0.1 && \
echo "###  step2:finished DMR analysis"  && rm -f "$calling_dir/$name.DMR.raw.txt" "$calling_dir/$name.DMR.log" &

metilene -a g1 -b g2 -f 3 -t $threads $calling_dir/$name.inputf.txt > "$calling_dir/$name.DMC.raw.txt" 2>"$calling_dir/$name.DMC.log" && \
sort --temporary-directory="./" --parallel=$threads -k1,1V -k2,2n "$calling_dir/$name.DMC.raw.txt" > "$calling_dir/$name.DMC.bed" && \
less $calling_dir/$name.DMC.bed |awk '($5>0.1 || $5<-0.1) && $7<0.05'|awk '{print $1,$2,$7,$5,$6,$9,$10}' OFS='\t' > $calling_dir/$name.DMC_pval.0.05.out && \
echo "###  step2:finished DMC analysis" && rm -f "$calling_dir/$name.DMC.raw.txt" "$calling_dir/$name.DMC.log" &

wait

end_time=$(date +%s)
echo "end time: $(date '+%H:%M:%S')"
duration=$((end_time - start_time))
hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))
seconds=$((duration % 60))
printf "all processes time: %02d:%02d:%02d\n" $hours $minutes $seconds

echo "##############################   finished  ############################################"
