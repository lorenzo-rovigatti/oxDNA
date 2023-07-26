#!/bin/bash

Nbp=17
Nb=$(($Nbp*2))

echo "{" >> op.txt
echo "order_parameter = bond" >> op.txt
echo "name = all_native_bonds" >> op.txt

for ((i=0; i< $Nbp; i++)); do
 n=$(($i+1))
 paired=$(($Nb-1-$i))
 echo "pair${n} = $i, ${paired}" >> op.txt
done

echo "}" >> op.txt
