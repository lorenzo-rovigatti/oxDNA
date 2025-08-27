#!/bin/bash -l

path=$(pwd)
oxDNA_path=/users/yqb22156/oxdna

nreps=1
box_size=20
NTs=8

seqs=("ATTAAATT" "GAATTTAA" "TAACCTAT" "AGTCTCTT" "CGATGCTA" "TGGTCGTT" "GTGGCAGC" "CGCAGCGT" "GCGCGCGA" "GCGCGCGC" )
mm_id=( 0 1 2 0 1 2 0 1 )
OP_files=( "op_n8.txt" "op_n8.txt" "op_n8.txt" "op_n8.txt" "op_n8.txt" "op_n8.txt" "op_n8.txt" "op_n8.txt" "op_n8.txt" "op_n8.txt" )
W_files=( "wfile_n8.txt" "wfile_n8.txt" "wfile_n8.txt" "wfile_n8.txt" "wfile_n8.txt" "wfile_n8.txt" "wfile_n8.txt" "wfile_n8.txt" "wfile_n8.txt" "wfile_n8.txt"  )
Nseqs=${#seqs[@]}
nbps=( 8 8 8 8 8 8 8 8 8 8 )
mTs=( 22 30 37 39 41 49 49 56 65 66 )

for ((l=0; l < ${Nseqs}; l++)); do
	for ((m=0; m < ${nreps}; m++)); do

		rep_path=${path}/Seq${l}/Rep${m}

		mkdir -p ${rep_path}

		cp ${path}/input_PT_VMMC ${rep_path}
        sub_script=PT_VMMC_seq${l}_rep${m}_n8.sh
        cp ${path}/PT_VMMC_seq0_rep0_n8.sh ${rep_path}/${sub_script}
		cp ${path}/${OP_files[$l]} ${rep_path}
		cp ${path}/${W_files[$l]} ${rep_path}
		cp ${path}/oxDNA_sequence_dependent_parameters_in.txt ${rep_path}
		cd ${rep_path}
		echo "double ${seqs[l]}" > gen.txt
		#python3 ${oxDNA_path}/utils/generate-sa.py ${box_size} gen.txt
		bash ../generate_tmm.sh ${box_size} ${mm_id[$l]}
        for ((k=0; k < ${NTs}; k++)); do
            cp generated.dat generated.dat${k}
            cp ${OP_files[$l]} ${OP_files[$l]}$k
            #cp ${W_files[$l]} ${W_files[$l]}$k
        done
		sed -i "s|seed = .*|seed = ${RANDOM}|g" input_PT_VMMC
		sed -i "s|T = .*|T = ${mTs[$l]}C|g" input_PT_VMMC

        sed -i "s|l=.*|l=${l}|g" ${sub_script}
        sed -i "s|m=.*|m=${m}|g" ${sub_script}
        sed -i "s|path=.*|path=${path}|g" ${sub_script}
        sed -i "s|job-name=.*|job-name=PTVMMC_s${l}r${m}_n8|g" ${sub_script}

		extr_Ts="extrapolate_hist ="
        sim_Ts="pt_temp_list ="
        in_n=$(( (${mTs[$l]}-2)/4 ))
        if [ $in_n -gt 3 ]; then
            in_n=3
        fi


        counts=-1
        fin_n=$(( $NTs-1-${in_n} ))
        for ((n=-${in_n}; n<=${fin_n}; n++)); do
            counts=$(($counts+1))
            ex_T=$(( $n*4 + ${mTs[$l]} -2 ))
            tmp=$(( $n*4 -2 ))
            delta=""
            #echo "tmp = ${tmp}"
            if [ ${tmp} -lt 0 ]; then
                tmp=$((-${tmp}))
                delta="m${tmp}"
            else
                delta="p${tmp}"
            fi
            #echo "delta = ${delta}"
            cp ${path}/Weights_n8/wfile_n8_${delta}.txt ${rep_path}/wfile_n8.txt${counts}
            if [ ${ex_T} = 0 ]; then
                ex_T=$(( ${ex_T}+1 ))
            fi
            if [ $n -lt ${fin_n} ]; then
                sim_Ts="${sim_Ts} ${ex_T}.0C,"
            else
                sim_Ts="${sim_Ts} ${ex_T}.0C"
            fi
        done

        in_n=$(( $in_n*2 ))
        fin_n=$(( $fin_n*2 ))
        for ((n=-${in_n}; n<=${fin_n}; n++)); do
            ex_T=$(( $n*2 + ${mTs[$l]} -1 ))
            if [ ${ex_T} = 0 ]; then
                ex_T=$(( ${ex_T}+1 ))
            fi
            if [ $n -lt ${fin_n} ]; then
                extr_Ts="${extr_Ts} ${ex_T}.0C,"
            else
                extr_Ts="${extr_Ts} ${ex_T}.0C"
            fi
        done


		sed -i "s|extrapolate_hist = .*|${extr_Ts}|g" input_PT_VMMC
        sed -i "s|pt_temp_list = .*|${sim_Ts}|g" input_PT_VMMC
		#${oxDNA_path}/build/bin/oxDNA input_VMMC_relax > out &

	done
done

wait
