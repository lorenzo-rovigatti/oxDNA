path=$(pwd)
nseq=10
nrep=1

for (( l=0; l<${nseq}; l++ )); do
    for (( m=0; m<${nrep}; m++ )); do
        cd ${path}/Seq${l}/Rep${m}/
        sbatch PT_VMMC_seq${l}_rep${m}_n8.sh
    done
done
