#!/usr/bin/bash

# List of targets
TARGETS=( \
    test/DNA/DSDNA8/MC \
    test/DNA/DSDNA8/MD \
    test/DNA/DSDNA8/VMMC \
    test/DNA/HAIRPIN18/VMMC_DOUBLE \
    test/DNA/HAIRPIN18/VMMC_DOUBLE \
    test/DNA/SSDNA15/MC \
    test/DNA/SSDNA15/MD \
    test/DNA/SSDNA15/VMMC \
    test/LJ \
    test/THERMOSTATS/BUSSI \
    test/THERMOSTATS/JOHN \
    test/THERMOSTATS/LANGEVIN \
)

# Find oxDNA and link to temporary directory
ln $(find . -name oxDNA) /tmp/oxDNA
OXDNA_BIN=/tmp/oxDNA
echo "OXDNA_BIN: ${OXDNA_BIN}"

# Run all of the test scripts
for target in ${TARGETS[@]}; do
pushd $target >& /dev/null
echo "RUNNING: ${target}"
$OXDNA_BIN quick_input
popd >& /dev/null
echo "FINISHED!"
done

# remove the linker
rm /tmp/oxDNA
