if [[ $# -ne 1 ]]; then
    echo "Error: Exactly 1 argument required" >&2
    echo "Usage: $0 box_size" >&2
    exit 1
fi

box_size=$1

path=$(pwd)

oxdna_path=/home/users/yqb22156/oxdna

python3 ${oxdna_path}/utils/generate-sa.py ${box_size} seq.txt

#topology
awk 'NR==1{$1=$1-1} {print}' generated.top > output.txt && mv output.txt generated.top #subtract 1 from first entry (number of nucleotides)
sed -i '$d' generated.top #remove last line
awk 'NR==FNR{if(NR>1)print prev; prev=$0; next} END{$NF=-1; print}' generated.top generated.top > tmp && mv tmp generated.top #change last entry to -1

#coordinates
sed -i '$d' generated.dat #remove last line
