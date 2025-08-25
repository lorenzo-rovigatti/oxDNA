if [[ $# -ne 2 ]]; then
    echo "Error: Exactly 2 arguments required" >&2
    echo "Usage: $0 box_size mismatch_id(=0,1,2)" >&2
    exit 1
fi

box_size=$1
mm_id=$2

path=$(pwd)

oxdna_path=/home/users/yqb22156/oxdna

python3 ${oxdna_path}/utils/generate-sa.py $box_size seq.txt

#topology
mm_type=$(tail -n1 generated.top | awk '{print $2}')

declare -a c_types

if [ ${mm_type} = 'A' ]; then
	c_types=('G' 'C' 'T')
elif [ ${mm_type} = 'G' ]; then
	c_types=('A' 'C' 'T')
elif [ ${mm_type} = 'C' ]; then
	c_types=('A' 'G' 'T')
else
	c_types=('A' 'G' 'C')
fi


awk -v new="${c_types[${mm_id}]}" '
NR>1 {print prev}
{prev=$0}
END {
    $2=new
    print
}' generated.top > temp.txt && mv temp.txt generated.top

