#!/usr/bin/bash

# Find oxDNA and link to temporary directory
ln $(find . -name oxDNA) /tmp/oxDNA
OXDNA_BIN=/tmp/oxDNA
echo "OXDNA_BIN: ${OXDNA_BIN}"

# Run one test script to check that the program compiles
pushd test/DNA/FAST >& /dev/null
$OXDNA_BIN quick_input
popd >& /dev/null

# remove the linker
rm /tmp/oxDNA
