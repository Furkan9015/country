#!/bin/bash

set -x

# Load environment
source ~/.bashrc
source ~/upmem/upmem_env.sh

# Activate openjdl for Java but use system Python
eval "$(micromamba shell hook --shell bash)"
micromamba activate openjdl

# Make sure we use system python and java from openjdl
export PATH="/usr/bin:$PATH"

# Verify tools are available
echo "Python: $(which python)"
echo "Java: $(which java)"
python --version
java -version

GENOMEE=$(realpath $1)
CHRALL=$(realpath $2)
RESULT=$(pwd)/result.txt

rm -rf ${RESULT}

for ((c = 1; c < 25; c++))
do
    echo $c >> $RESULT
    date >> $RESULT
    rm -rf chr${c}_640
    mkdir -p chr${c}_640
    pushd chr${c}_640
    cp ${CHRALL}/chr${c}.fasta .
    cp ${CHRALL}/sel_var.py .
    cp ${CHRALL}/GenomeParser.py .
    cp ${CHRALL}/combineFastq.py .
    cat ${CHRALL}/chrallvars.vcf | grep -E "^${c}\>" > chr${c}vars.vcf
    date >> $RESULT
    time ${GENOMEE}/tests/simchr_mono.bash ${c}
    date >> $RESULT
    time ${GENOMEE}/build/host/upvc -i chr${c} -g index -n 640
    date >> $RESULT
    time ${GENOMEE}/build/host/upvc -i chr${c} -g map
    date >> $RESULT
    python3 ${GENOMEE}/tests/compareVCF.py *.vcf >> ${RESULT}
    popd
done
