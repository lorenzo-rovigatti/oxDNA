#!/bin/bash

echo "Testing align..."
if
    python ../src/oxDNA_analysis_tools/align.py minitraj.dat aligned.dat 2>&1 >/dev/null | grep -y "ERROR"
then
    echo "AN ERROR OCCURED"
else
    echo "OK"
fi

echo ""

echo "Testing backbone_flexibility..."
if
    python ../src/oxDNA_analysis_tools/backbone_flexibility.py rna_tile.top minitraj.dat 2>&1 >/dev/null | grep -y "ERROR"
then
    echo "AN ERROR OCCURED"
else
    echo "OK"
fi

echo ""

echo "Testing bond_analysis..."
if
    python ../src/oxDNA_analysis_tools/bond_analysis.py -p2 input_rna minitraj.dat pairs.txt 2>&1 >/dev/null | grep -y "ERROR"
then
    echo "AN ERROR OCCURED"
else
    echo "OK"
fi

echo ""

echo "Testing mean and deviations..."
if
    python  ../src/oxDNA_analysis_tools/mean.py -p 2 -d devs.json minitraj.dat 2>&1 >/dev/null | grep -y "ERROR"
then
    echo "AN ERROR OCCURED"
else
    echo "OK"
fi

echo ""

echo "Testing centroid with indexing..."
if
    python ../src/oxDNA_analysis_tools/centroid.py -i index.txt mean.dat minitraj.dat 2>&1 >/dev/null | grep -y "ERROR"
then
    echo "AN ERROR OCCURED"
else
    echo "OK"
fi

echo ""

echo "Testing contact_map..."
if
    python ../src/oxDNA_analysis_tools/contact_map.py minitraj.dat 2>&1 >/dev/null | grep -y "ERROR"
then
    echo "AN ERROR OCCURED"
else
    echo "OK"
fi

echo ""

echo "Testing distance and clustering (this one takes a while because of the plot)..."
if
    python ../src/oxDNA_analysis_tools/distance.py -c -i minitraj.dat 1 3 5 67 34 56 2>&1 >/dev/null | grep -y "ERROR"
then
    echo "AN ERROR OCCURED"
else
    echo "OK"
fi

echo ""

echo "Testing duplex_finder..."
if
    python ../src/oxDNA_analysis_tools/duplex_finder.py input_rna minitraj.dat 2>&1 >/dev/null | grep -y "ERROR"
then
    echo "AN ERROR OCCURED"
else
    echo "OK"
fi

echo ""

echo "Testing duplex_angle_plotter..."
if
    python ../src/oxDNA_analysis_tools/duplex_angle_plotter.py -o angle.png -i angles.txt 7 37 2>&1 >/dev/null | grep -y "ERROR"
then
    echo "AN ERROR OCCURED"
else
    echo "OK"
fi

echo ""

echo "Testing generate_force..."
if
    python ../src/oxDNA_analysis_tools/generate_force.py -f gen_pairs.txt input_rna minitraj.dat 2>&1 >/dev/null | grep -y "ERROR"
then
    echo "AN ERROR OCCURED"
else
    echo "OK"
fi

echo ""

echo "Testing minify..."
if
    python ../src/oxDNA_analysis_tools/minify.py -a -d2 minitraj.dat min.dat 2>&1 >/dev/null | grep -y "ERROR"
then
    echo "AN ERROR OCCURED"
else
    echo "OK"
fi

echo ""

echo "Testing multidimensional_scaling_mean..."
if
    python ../src/oxDNA_analysis_tools/multidimensional_scaling_mean.py minitraj.dat -o meanM -d devsM 2>&1 >/dev/null | grep -y "ERROR"
then
    echo "AN ERROR OCCURED"
else
    echo "OK"
fi

echo ""

echo "Testing output_bonds..."
if
    python ../src/oxDNA_analysis_tools/output_bonds.py -v energy.json input_rna minitraj.dat 2>&1 >/dev/null | grep -y "ERROR"
then
    echo "AN ERROR OCCURED"
else
    echo "OK"
fi

echo ""

echo "Testing pca"
if
    python ../src/oxDNA_analysis_tools/pca.py minitraj.dat mean.dat pca.json 2>&1 >/dev/null | grep -y "ERROR"
then
    echo "AN ERROR OCCURED"
else
    echo "OK"
fi

echo "Testing subset_trajectory..."
if
    python ../src/oxDNA_analysis_tools/subset_trajectory.py minitraj.dat rna_tile.top -i index.txt sub.dat 2>&1 >/dev/null | grep -y "ERROR"
then
    echo "AN ERROR OCCURED"
else
    echo "OK"
fi

echo ""

echo "Testing superimpose..."
if
    python ../src/oxDNA_analysis_tools/superimpose.py mean.dat centroid.dat 2>&1 >/dev/null | grep -y "ERROR"
then
    echo "AN ERROR OCCURED"
else
    echo "OK"
fi

echo ""
