#!/bin/bash

echo "Testing oat"
echo "-----------"

if [ $# -gt 0 ];
then
    PYTHON=$1
else
    PYTHON="python"
fi

test() {
    tmpfile=$(mktemp)
    $PYTHON $1 &> $tmpfile
    local code_res=$?
    grep_out=$(grep -i "ERROR" $tmpfile)
    local grep_res=$?

    if [ $code_res -eq 0 ] && [ $grep_res -eq 1 ]
    then
        echo " OK"
    else
        echo ""
        echo $grep_out
        echo "--> AN ERROR OCCURRED <--"
    fi
}

echo -n "Testing align..."
test "../src/oxDNA_analysis_tools/align.py minitraj.dat aligned.dat"

echo -n "Testing backbone_flexibility..."
test "../src/oxDNA_analysis_tools/backbone_flexibility.py rna_tile.top minitraj.dat"

echo -n "Testing bond_analysis..."
test "../src/oxDNA_analysis_tools/bond_analysis.py -p2 input_rna minitraj.dat pairs.txt"

echo -n "Testing mean and deviations..."
test "../src/oxDNA_analysis_tools/mean.py -p 2 -d devs.json minitraj.dat"

echo -n "Testing centroid with indexing..."
test "../src/oxDNA_analysis_tools/centroid.py -i index.txt mean.dat minitraj.dat"

echo -n "Testing contact_map..."
test "../src/oxDNA_analysis_tools/contact_map.py minitraj.dat"

echo -n "Testing distance and clustering (this one takes a while because of the plot)..."
test "../src/oxDNA_analysis_tools/distance.py -c -i minitraj.dat 1 3 5 67 34 56"

echo -n "Testing duplex_finder..."
test "../src/oxDNA_analysis_tools/duplex_finder.py input_rna minitraj.dat"

echo -n "Testing duplex_angle_plotter..."
test "../src/oxDNA_analysis_tools/duplex_angle_plotter.py -o angle.png -i angles.txt 7 37"

echo -n "Testing generate_force..."
test "../src/oxDNA_analysis_tools/generate_force.py -f gen_pairs.txt input_rna minitraj.dat"

echo -n "Testing minify..."
test "../src/oxDNA_analysis_tools/minify.py -a -d2 minitraj.dat min.dat"

echo -n "Testing multidimensional_scaling_mean..."
test "../src/oxDNA_analysis_tools/multidimensional_scaling_mean.py minitraj.dat -o meanM -d devsM"

echo -n "Testing output_bonds..."
test "../src/oxDNA_analysis_tools/output_bonds.py -v energy.json input_rna minitraj.dat"

echo -n "Testing pca"
test "../src/oxDNA_analysis_tools/pca.py minitraj.dat mean.dat pca.json"

echo -n "Testing subset_trajectory..."
test "../src/oxDNA_analysis_tools/subset_trajectory.py minitraj.dat rna_tile.top -i index.txt sub.dat"

echo -n "Testing superimpose..."
test "../src/oxDNA_analysis_tools/superimpose.py mean.dat centroid.dat"

echo ""