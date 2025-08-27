if [[ $# -ne 1 ]]; then
    echo "Error: Exactly 1 argument required" >&2
    echo "Usage: $0 mm_pos" >&2
    exit 1
fi

pos=$1

echo "Will generate mismatch after the first $1 bps"

path=$(pwd)

oxdna_path=/home/users/yqb22156/oxdna

python3 ${oxdna_path}/utils/generate-sa.py 20 seq.txt


line_num=$(($pos+2))

#echo ${line_num}

#topology
mm_type=$(awk -v line="$line_num" 'NR==line {print $2}' generated.top)

#echo $mm_type

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

#echo ${c_types[@]}


awk -v new="${c_types[0]}" -v pos="$line_num" 'NR==pos{$2=new}1' OFS=' ' generated.top > temp && mv temp generated_mm0.top
awk -v new="${c_types[1]}" -v pos="$line_num" 'NR==pos{$2=new}1' OFS=' ' generated.top > temp && mv temp generated_mm1.top
awk -v new="${c_types[2]}" -v pos="$line_num" 'NR==pos{$2=new}1' OFS=' ' generated.top > temp && mv temp generated_mm2.top
