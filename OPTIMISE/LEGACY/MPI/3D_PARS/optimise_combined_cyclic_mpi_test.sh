main_path_cycle=$(pwd)
opti_path_cycle=/home/users/yqb22156/Desktop/oxDNA/Parametrisation/Minimisation_ox3/OPTIMISE/MPI
echo "Usage: "
echo "bash optimise.sh config_file_melting c_geo1 c_geo2 c_geo3 c_geo4"

if [ "$#" -ne 5 ]; then
	echo "Illegal number of parameters"
	exit 1
fi

config_m=$1 #configuration file
config_g0=$2
config_g1=$3
config_g2=$4
config_g3=$5

mkdir ${main_path_cycle}/Step0_melting
cp ${main_path_cycle}/$config_m Step0_melting/
cp ${main_path_cycle}/op_n8* Step0_melting/
cp ${main_path_cycle}/wfile* $main_path_cycle/Step0_melting/
cp ${main_path_cycle}/wfile* $main_path_cycle/Step0_melting/

cp ${main_path_cycle}/oxDNA_sequence_dependent_parameters.txt ${main_path_cycle}/Step0_melting
cp ${main_path_cycle}/oxDNA_sequence_dependent_parameters_in.txt ${main_path_cycle}/Step0_melting


cd $main_path_cycle/Step0_melting
bash $opti_path_cycle/optimise_multi_melting_mpi.sh $config_m > MB_OUT

wait
###############################################
cd $main_path_cycle/

mkdir -p Step1_geo/0/
mkdir -p Step1_geo/1/
mkdir -p Step1_geo/2/
mkdir -p Step1_geo/3/

par_file=$main_path_cycle/Step0_melting/Step1/oxDNA_sequence_dependent_parameters_fin.txt
cp $par_file $main_path_cycle/Step1_geo/oxDNA_sequence_dependent_parameters_in.txt

tail -n 20 $par_file > $main_path_cycle/Step1_geo/oxDNA_sequence_dependent_parameters.txt
awk '/STCK_FACT_EPS|FENE_DELTA_|HYDR_THETA4_|CRST_K_|CRST_K_/ {print $0}' $par_file >> $main_path_cycle/Step1_geo/oxDNA_sequence_dependent_parameters.txt


cp $main_path_cycle/$config_g0 $main_path_cycle/Step1_geo/0/
cp $main_path_cycle/Step1_geo/oxDNA_sequence_dependent_parameters.txt $main_path_cycle/Step1_geo/0/
cp $main_path_cycle/Step1_geo/oxDNA_sequence_dependent_parameters_in.txt $main_path_cycle/Step1_geo/0/

cd $main_path_cycle/Step1_geo/0/
bash $opti_path_cycle/optimise_multi_mpi.sh $config_g0 > GB_OUT &

cp $main_path_cycle/$config_g1 $main_path_cycle/Step1_geo/1/
cp $main_path_cycle/Step1_geo/oxDNA_sequence_dependent_parameters.txt $main_path_cycle/Step1_geo/1/
cp $main_path_cycle/Step1_geo/oxDNA_sequence_dependent_parameters_in.txt $main_path_cycle/Step1_geo/1/

cd $main_path_cycle/Step1_geo/1/
bash $opti_path_cycle/optimise_multi_mpi.sh $config_g1 > GB_OUT &

cp $main_path_cycle/$config_g2 $main_path_cycle/Step1_geo/2/
cp $main_path_cycle/Step1_geo/oxDNA_sequence_dependent_parameters.txt $main_path_cycle/Step1_geo/2/
cp $main_path_cycle/Step1_geo/oxDNA_sequence_dependent_parameters_in.txt $main_path_cycle/Step1_geo/2/

cd $main_path_cycle/Step1_geo/2/
bash $opti_path_cycle/optimise_multi_mpi.sh $config_g2 > GB_OUT &

cp $main_path_cycle/$config_g3 $main_path_cycle/Step1_geo/3/
cp $main_path_cycle/Step1_geo/oxDNA_sequence_dependent_parameters.txt $main_path_cycle/Step1_geo/3/
cp $main_path_cycle/Step1_geo/oxDNA_sequence_dependent_parameters_in.txt $main_path_cycle/Step1_geo/3/ 

cd $main_path_cycle/Step1_geo/3/
bash $opti_path_cycle/optimise_multi_mpi.sh $config_g3 > GB_OUT &

wait
#################################
cd $main_path_cycle/

mkdir Step2_melting

cp $main_path_cycle/$config_m Step2_melting/
cp $main_path_cycle/op_n8* Step2_melting/
cp $main_path_cycle/Step0_melting/Step1/wfile* $main_path_cycle/Step2_melting/

cd Step2_melting

cp $main_path_cycle/Step1_geo/0/oxDNA_sequence_dependent_parameters_fin.txt $main_path_cycle/Step2_melting/oxDNA_sequence_dependent_parameters_in.txt
awk 'NR>101' $main_path_cycle/Step1_geo/1/oxDNA_sequence_dependent_parameters_fin.txt >>  $main_path_cycle/Step2_melting/oxDNA_sequence_dependent_parameters_in.txt
awk 'NR>101' $main_path_cycle/Step1_geo/2/oxDNA_sequence_dependent_parameters_fin.txt >>  $main_path_cycle/Step2_melting/oxDNA_sequence_dependent_parameters_in.txt
awk 'NR>101' $main_path_cycle/Step1_geo/3/oxDNA_sequence_dependent_parameters_fin.txt >>  $main_path_cycle/Step2_melting/oxDNA_sequence_dependent_parameters_in.txt

awk 'NR>20' $main_path_cycle/Step2_melting/oxDNA_sequence_dependent_parameters_in.txt >>  $main_path_cycle/Step2_melting/oxDNA_sequence_dependent_parameters.txt

bash $opti_path_cycle/optimise_multi_melting_mpi.sh $config_m > MB_OUT

wait
####################

cd $main_path_cycle/

mkdir -p Step3_geo/0/
mkdir -p Step3_geo/1/
mkdir -p Step3_geo/2/
mkdir -p Step3_geo/3/

par_file=$main_path_cycle/Step2_melting/Step1/oxDNA_sequence_dependent_parameters_fin.txt
cp $par_file $main_path_cycle/Step3_geo/oxDNA_sequence_dependent_parameters_in.txt

tail -n 20 $par_file > $main_path_cycle/Step3_geo/oxDNA_sequence_dependent_parameters.txt
awk '/STCK_FACT_EPS|FENE_DELTA_|HYDR_THETA4_|CRST_K_|CRST_K_/ {print $0}' $par_file >> $main_path_cycle/Step3_geo/oxDNA_sequence_dependent_parameters.txt

cp $main_path_cycle/$config_g0 $main_path_cycle/Step3_geo/0/
cp $main_path_cycle/Step3_geo/oxDNA_sequence_dependent_parameters.txt $main_path_cycle/Step3_geo/0/
cp $main_path_cycle/Step3_geo/oxDNA_sequence_dependent_parameters_in.txt $main_path_cycle/Step3_geo/0/

cd $main_path_cycle/Step3_geo/0/
bash $opti_path_cycle/optimise_multi_mpi.sh $config_g0 > GB_OUT &

cp $main_path_cycle/$config_g1 $main_path_cycle/Step3_geo/1/
cp $main_path_cycle/Step3_geo/oxDNA_sequence_dependent_parameters.txt $main_path_cycle/Step3_geo/1/
cp $main_path_cycle/Step3_geo/oxDNA_sequence_dependent_parameters_in.txt $main_path_cycle/Step3_geo/1/

cd $main_path_cycle/Step3_geo/1/
bash $opti_path_cycle/optimise_multi_mpi.sh $config_g1 > GB_OUT &

cp $main_path_cycle/$config_g2 $main_path_cycle/Step3_geo/2/
cp $main_path_cycle/Step3_geo/oxDNA_sequence_dependent_parameters.txt $main_path_cycle/Step3_geo/2/
cp $main_path_cycle/Step3_geo/oxDNA_sequence_dependent_parameters_in.txt $main_path_cycle/Step3_geo/2/

cd $main_path_cycle/Step3_geo/2/
bash $opti_path_cycle/optimise_multi_mpi.sh $config_g2 > GB_OUT &

cp $main_path_cycle/$config_g3 $main_path_cycle/Step3_geo/3/
cp $main_path_cycle/Step3_geo/oxDNA_sequence_dependent_parameters.txt $main_path_cycle/Step3_geo/3/
cp $main_path_cycle/Step3_geo/oxDNA_sequence_dependent_parameters_in.txt $main_path_cycle/Step3_geo/3/ 

cd $main_path_cycle/Step3_geo/3/
bash $opti_path_cycle/optimise_multi_mpi.sh $config_g3 > GB_OUT &

wait

##################################
