srcdir=/XYFS01/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_2/BIGDATA2/gzfezx_shhli_2/USER/luozhenyu/script/rrbs
python $srcdir/tools/createpipeline.py \
        -i undo_89_smaples.info \
        -c $srcdir/config.yaml \
        -r $srcdir/run_snakemake.sh \
        -j 8 \
        -p ./
