#!/bin/bash
# run_benchmark.sh
#
# Benchmarks the new direct CUDAEdgeList (CUDA_list = edge) against the
# existing Verlet list (CUDA_list = verlet, with and without use_edge) on
# GPU, and Verlet vs Cells on CPU for reference.
# Uses DNA2 systems from ErikPoppleton/oxDNA_performance.
# Records wall time, steps/s, and peak VRAM usage per run.
#
# Usage:
#   ./run_benchmark.sh [build_dir] [steps] [bench_data_dir]
#
# Arguments:
#   build_dir       path to the oxDNA build directory (default: ../build)
#   steps           number of MD steps (default: 100000)
#   bench_data_dir  path to a checkout of ErikPoppleton/oxDNA_performance
#                   (default: ./oxDNA_performance next to this script)

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${1:-${SCRIPT_DIR}/../build}"
STEPS="${2:-100000}"
BENCH_DATA="${3:-${SCRIPT_DIR}/oxDNA_performance}"
RESULTS="${SCRIPT_DIR}/results"
OXDNA="${BUILD_DIR}/bin/oxDNA"

mkdir -p "${RESULTS}"

if [ ! -x "${OXDNA}" ]; then
    echo "ERROR: oxDNA binary not found at ${OXDNA}"
    echo "Build the project first (see docs/source/install.md), e.g.:"
    echo "  cmake -S . -B build -DCUDA=ON && cmake --build build -j"
    exit 1
fi

SIZES="N64 N512 N4096"

GPU_AVAILABLE=0
if command -v nvidia-smi &>/dev/null && nvidia-smi &>/dev/null; then
    GPU_AVAILABLE=1
    echo "GPU detected: $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
fi

echo "case,time_ms,steps_per_sec,vram_mb" > "${RESULTS}/summary.csv"

# _poll_vram PEAK_FILE
# Background VRAM poller — writes running max to PEAK_FILE every 0.25 s.
_poll_vram() {
    local peak_file="$1"
    while true; do
        local cur
        cur=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ' ')
        local prev
        prev=$(cat "${peak_file}" 2>/dev/null || echo 0)
        if [ -n "${cur}" ] && [ "${cur}" -gt "${prev}" ] 2>/dev/null; then
            echo "${cur}" > "${peak_file}"
        fi
        sleep 0.25
    done
}

# run_case NAME INPUT_FILE OUT_FILE LOG_DIR IS_GPU
run_case() {
    local name="$1"
    local input_file="$2"
    local out_file="$3"
    local log_dir="$4"
    local is_gpu="${5:-0}"

    echo -n "  ${name} ... "

    # VRAM: capture baseline and start poller for GPU runs
    local vram_str="N/A"
    local peak_file poller_pid baseline_vram
    if [ "${is_gpu}" -eq 1 ] && [ "${GPU_AVAILABLE}" -eq 1 ]; then
        baseline_vram=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ' ')
        peak_file=$(mktemp)
        echo "${baseline_vram}" > "${peak_file}"
        _poll_vram "${peak_file}" &
        poller_pid=$!
    fi

    local t_start t_end elapsed ret
    t_start=$(date +%s%N)
    "${OXDNA}" "${input_file}" > "${out_file}" 2>&1
    ret=$?
    t_end=$(date +%s%N)
    elapsed=$(( (t_end - t_start) / 1000000 ))

    # Stop VRAM poller and compute delta
    if [ "${is_gpu}" -eq 1 ] && [ "${GPU_AVAILABLE}" -eq 1 ]; then
        kill "${poller_pid}" 2>/dev/null
        wait "${poller_pid}" 2>/dev/null || true
        local peak_vram
        peak_vram=$(cat "${peak_file}" 2>/dev/null || echo "${baseline_vram}")
        rm -f "${peak_file}"
        vram_str=$(( peak_vram - baseline_vram ))
    fi

    if [ "${ret}" -ne 0 ]; then
        echo "FAILED"
        tail -5 "${out_file}"
        return
    fi

    # parse SimBackend total time from the oxDNA log file
    local wall_sec sps
    wall_sec=$(grep "> SimBackend" "${log_dir}/log.dat" 2>/dev/null \
               | grep -oP '^\s*> SimBackend\s+\K[\d.]+' | head -1 || true)
    if [ -n "${wall_sec}" ] && [ "${wall_sec}" != "0" ]; then
        sps=$(awk "BEGIN { printf \"%.0f\", ${STEPS} / ${wall_sec} }")
    else
        sps="N/A"
    fi

    echo "${elapsed}ms (wall) / sim=${wall_sec:-?}s — ${sps} steps/s  VRAM=${vram_str} MB"
    echo "${name},${elapsed},${sps},${vram_str}" >> "${RESULTS}/summary.csv"
}

# make_input OUTPUT_FILE EXTRA_PARAMS TOP CONF TMPDIR
make_input() {
    local output="$1"
    local extra="$2"
    local top="$3"
    local conf="$4"
    local tmpd="$5"

    mkdir -p "${tmpd}"
    cat > "${output}" <<EOFINPUT
sim_type = MD
steps = ${STEPS}
dt = 0.003
ensemble = nvt
T = 20C
salt_concentration = 1.0
thermostat = brownian
newtonian_steps = 103
diff_coeff = 2.5
interaction_type = DNA2
use_average_seq = 1
verlet_skin = 0.5
external_forces = 0
topology = ${top}
conf_file = ${conf}
lastconf_file = ${tmpd}/last.dat
trajectory_file = ${tmpd}/traj.dat
energy_file = ${tmpd}/energy.dat
log_file = ${tmpd}/log.dat
print_conf_interval = ${STEPS}
print_energy_every = ${STEPS}
time_scale = linear
refresh_vel = 1
restart_step_counter = 1
no_stdout_energy = 0
${extra}
EOFINPUT
}

for SIZE in ${SIZES}; do
    BENCH="${BENCH_DATA}/${SIZE}"
    if [ ! -d "${BENCH}" ]; then
        echo "Skipping ${SIZE}: directory not found"
        continue
    fi

    TOP=$(ls "${BENCH}"/topology_*.top 2>/dev/null | head -1 || true)
    CONF=$(ls "${BENCH}"/init_conf_*.dat 2>/dev/null | head -1 || true)
    if [ -z "${CONF}" ]; then
        CONFZIP=$(ls "${BENCH}"/init_conf_*.zip 2>/dev/null | head -1 || true)
        if [ -n "${CONFZIP}" ]; then
            echo "  Unzipping ${CONFZIP}..."
            unzip -q "${CONFZIP}" -d "${BENCH}"
            CONF=$(ls "${BENCH}"/init_conf_*.dat 2>/dev/null | head -1 || true)
        fi
    fi
    if [ -z "${TOP}" ] || [ -z "${CONF}" ]; then
        echo "  ${SIZE}: missing topology or conf, skipping"
        continue
    fi

    echo ""
    echo "=== ${SIZE} ($(wc -l < "${TOP}") lines topology) ==="

    # ----------------------------------------------------------------
    # CPU benchmarks (no VRAM tracking): Verlet vs Cells, for reference.
    # Both are pre-existing oxDNA list types, unrelated to this PR's
    # CUDAEdgeList; included so the README graph has a CPU baseline.
    # ----------------------------------------------------------------
    echo "--- CPU ---"
    for LT in verlet cells; do
        TMPD="${RESULTS}/tmp_${SIZE}_cpu_${LT}"
        INP="${RESULTS}/${SIZE}_cpu_${LT}.in"
        OUT="${RESULTS}/${SIZE}_cpu_${LT}.log"
        make_input "${INP}" "backend = CPU
list_type = ${LT}" "${TOP}" "${CONF}" "${TMPD}"
        run_case "${SIZE}_cpu_${LT}" "${INP}" "${OUT}" "${TMPD}" 0
    done

    # ----------------------------------------------------------------
    # GPU benchmarks (with VRAM tracking)
    # ----------------------------------------------------------------
    if [ "${GPU_AVAILABLE}" -eq 1 ]; then
        echo "--- GPU ---"

        # baseline: plain Verlet list, per-particle force kernels
        TMPD="${RESULTS}/tmp_${SIZE}_gpu_verlet"
        INP="${RESULTS}/${SIZE}_gpu_verlet.in"
        OUT="${RESULTS}/${SIZE}_gpu_verlet.log"
        make_input "${INP}" "backend = CUDA
backend_precision = mixed
CUDA_list = verlet
use_edge = 0" "${TOP}" "${CONF}" "${TMPD}"
        run_case "${SIZE}_gpu_verlet" "${INP}" "${OUT}" "${TMPD}" 1

        # existing approach: Verlet list + use_edge (matrix -> compress)
        TMPD="${RESULTS}/tmp_${SIZE}_gpu_verlet_edge"
        INP="${RESULTS}/${SIZE}_gpu_verlet_edge.in"
        OUT="${RESULTS}/${SIZE}_gpu_verlet_edge.log"
        make_input "${INP}" "backend = CUDA
backend_precision = mixed
CUDA_list = verlet
use_edge = 1
edge_n_forces = 1" "${TOP}" "${CONF}" "${TMPD}"
        run_case "${SIZE}_gpu_verlet_edge" "${INP}" "${OUT}" "${TMPD}" 1

        # new CUDAEdgeList: direct two-pass edge build, no N x max_neigh matrix
        TMPD="${RESULTS}/tmp_${SIZE}_gpu_edge"
        INP="${RESULTS}/${SIZE}_gpu_edge.in"
        OUT="${RESULTS}/${SIZE}_gpu_edge.log"
        make_input "${INP}" "backend = CUDA
backend_precision = mixed
CUDA_list = edge
use_edge = 1
edge_n_forces = 1" "${TOP}" "${CONF}" "${TMPD}"
        run_case "${SIZE}_gpu_edge" "${INP}" "${OUT}" "${TMPD}" 1
    else
        echo "--- GPU not available ---"
    fi
done

echo ""
echo "=== Done ==="
echo "Results in: ${RESULTS}/summary.csv"
echo "Run plot_results.py to generate a comparison graph."
echo ""
column -t -s, "${RESULTS}/summary.csv" 2>/dev/null || cat "${RESULTS}/summary.csv"
