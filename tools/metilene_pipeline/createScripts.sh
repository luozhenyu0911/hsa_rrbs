#!/bin/bash
name="$1"
# Check if directory exists and prompt user to recreate or exit
[ -d "$name" ] && {
    read -p "Directory '$name' exists. Recreate? (y/n): " -n 1 -r
    echo
    [[ $REPLY =~ ^[Yy]$ ]] && rm -rf "$name" || exit 1
}

mkdir $name && cd $name && mkdir CG && \
cp /XYFS02/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_1/Projects/20260514_S256rrbs/bin/config.conf CG/config.conf && \
# cp /XYFS02/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_1/Projects/20260514_S256rrbs/bin/config.conf CHG/config.conf && \
# cp /XYFS02/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_1/Projects/20260514_S256rrbs/bin/config.conf CHH/config.conf && \
echo "Created $name directory and copied config.conf files to each subdirectory."

if [ $# -ge 2 ] && [ "$3" = "run" ]; then

    echo "next step, you want to analyze the ${2}th column and $name in the pheno file"
    configf="/XYFS02/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_1/Projects/20260514_S256rrbs/bin/config.conf"
    script_py="/XYFS02/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_1/Projects/20260514_S256rrbs/bin/modify_config.py"
    python $script_py -c $configf -o CG/config.conf -n $name -col $2 \
    cd CG && yhbatch -N 1 -n 64 -p deimos /XYFS02/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_1/Projects/20260514_S256rrbs/bin/s1_call_DMR_DMC.sh && cd .. && \
    # python $script_py -c $configf -o CHG/config.conf -n $name -col $2 && \
    # cd CHG && yhbatch -N 1 -n 64 -p deimos /XYFS02/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_1/Projects/20260514_S256rrbs/bin/s1_call_DMR_DMC.sh && cd .. && \
    # python $script_py -c $configf -o CHH/config.conf -n $name -col $2 && \
    # cd CHH && yhbatch -N 1 -n 64 -p deimos /XYFS02/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_1/Projects/20260514_S256rrbs/bin/s1_call_DMR_DMC.sh && cd ..
fi














