#!/bin/bash 

for dir in `ls -d *`
do
cd $dir
 cat ds/trajectory.dat.FIT > ../$dir.FIT
 cat ss/trajectory.dat.FIT | tail -n+3 >> ../$dir.FIT 

cd ..
done
