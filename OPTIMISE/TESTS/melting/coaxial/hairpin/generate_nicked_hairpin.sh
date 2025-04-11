if [[ $# -ne 1 ]]; then
    echo "Error: Exactly 1 argument required" >&2
    echo "Usage: $0 mm_pos" >&2
    exit 1
fi

pos=$1

echo "Will insert nick after the first $1 bs"

path=$(pwd)

oxdna_path=/home/users/yqb22156/oxdna

python3 ${oxdna_path}/utils/generate-sa.py 20 seq.txt

cp generated_hairpin.dat generated.dat #replace generated.dat with hairpin initial configuration

line_num1=$(($pos+1))
line_num2=$(($pos+2))

#echo ${line_num}

#topology
awk 'NR==1{$2=$2+1} {print}' generated.top > output.txt && mv output.txt generated.top #add 1 from second entry (number of strands)
awk -v pos="$line_num1" 'NR==pos{$4=-1}1' OFS=' ' generated.top > temp && mv temp generated.top
awk -v pos="$line_num2" 'NR==pos{$3=-1}1' OFS=' ' generated.top > temp && mv temp generated.top

