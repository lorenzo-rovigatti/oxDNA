#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
FFS_SCRIPTS="${SCRIPT_DIR}/../scripts"

if [ -d FLUX ]; then
    echo "--> Entering FLUX"
    cd FLUX
    python3 "$FFS_SCRIPTS/ffs_flux.py" ffs.toml
    cd ..
fi

for d in $(find . -maxdepth 1 -type d -name 'SHOOT_*' | sort -V); do
    echo "--> Entering ${d#./}"
    cd "$d"
    python3 "$FFS_SCRIPTS/ffs_shoot.py" ffs.toml
    cd ..
done
