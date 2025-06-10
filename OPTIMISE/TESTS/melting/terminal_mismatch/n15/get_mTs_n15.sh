path=$(pwd)

path_prog=${path}/..

nseq=10
nrep=1
nreplicas=8

mv ${path}/measured_mTs_n15.txt ${path}/measured_mTs_n15_previus.txt

for ((i=0; i<${nseq}; i++)); do
    for ((j=0; j<${nrep}; j++)); do
        cd ${path}/Seq${i}/Rep${j}/
        python3 ${path_prog}/combine_histos.py ${nreplicas}
        python3 ${path_prog}/estimate_Tm_py3.py last_hist_comb.dat >> ${path}/measured_mTs_n15.txt
    done
done
