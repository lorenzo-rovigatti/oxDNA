#!/usr/bin/bash
OXDNA_BIN=./build/bin/oxDNA

$OXDNA_BIN test/DNA/DSDNA8/MC/quick_input
$OXDNA_BIN test/DNA/DSDNA8/MD/quick_input
$OXDNA_BIN test/DNA/DSDNA8/VMMC/quick_input
$OXDNA_BIN test/DNA/HAIRPIN18/VMMC_DOUBLE/quick_input
$OXDNA_BIN test/DNA/HAIRPIN18/VMMC_DOUBLE/quick_input
$OXDNA_BIN test/DNA/SSDNA15/MC/quick_input
$OXDNA_BIN test/DNA/SSDNA15/MD/quick_input
$OXDNA_BIN test/DNA/SSDNA15/VMMC/quick_input
$OXDNA_BIN test/LJ/quick_input
$OXDNA_BIN test/THERMOSTATS/BUSSI/quick_input
$OXDNA_BIN test/THERMOSTATS/JOHN/quick_input
$OXDNA_BIN test/THERMOSTATS/LANGEVIN/quick_input
