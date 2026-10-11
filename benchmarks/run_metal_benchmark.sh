#!/bin/bash
# oxDNA Metal-backend benchmark.
#
# Uses the duplex-box systems from ErikPoppleton/oxDNA_performance
# (128 - 524288 nucleotides, DNA2). Reports wall time per MD step, obtained by
# running two step counts and subtracting so the fixed start-up cost cancels;
# best of two repeats.
#
# Usage:
#   git clone https://github.com/ErikPoppleton/oxDNA_performance
#   (cd oxDNA_performance/N32768 && unzip init_conf_N32768.zip)
#   PERF=/path/to/oxDNA_performance \
#   OX_METAL=build_metal/bin/oxDNA OX_METAL_OMP=build_metal_omp/bin/oxDNA OX_CPU=build_cpu/bin/oxDNA \
#   benchmarks/run_metal_benchmark.sh > benchmarks/results.csv
#
# then:  python3 benchmarks/plot_benchmark.py benchmarks/results.csv benchmarks/plots
set -u

PERF="${PERF:?set PERF to the oxDNA_performance checkout}"
OXM="${OX_METAL:?set OX_METAL to a -DMETAL=ON oxDNA binary}"
OXO="${OX_METAL_OMP:-$OXM}"          # -DMETAL=ON -DUSE_OPENMP=ON binary (for hardmixed); falls back to OXM
OXC="${OX_CPU:-}"                    # -DDOUBLE=ON CPU binary; leave empty to skip CPU
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

echo "system,nuc,config,steps_short,steps_long,wall_short_s,wall_long_s,ms_per_step,total_1e6_steps_h,final_U"

# name   nuc      short long  run-cpu?
SYSTEMS=(
  "N8      128     400   4400  1"
  "N64     1024    400   4400  1"
  "N512    8192    300   3300  1"
  "N4096   65536   100   1100  1"
  "N32768  524288  40    240   0"
)

mkinput() {  # dir  backend-block  steps  sys
  cat > "$1/input" <<EOF
sim_type = MD
interaction_type = DNA2
salt_concentration = 1.0
use_average_seq = 1
T = 20C
dt = 0.003
verlet_skin = 0.5
thermostat = brownian
newtonian_steps = 103
diff_coeff = 2.5
seed = 12345
steps = $3
refresh_vel = 1
restart_step_counter = 1
topology = topology_$4.top
conf_file = $4_conf.dat
trajectory_file = tj.dat
energy_file = e.dat
print_energy_every = $3
print_conf_interval = 100000000
time_scale = linear
max_density_multiplier = 15
max_backbone_force = 10
$2
EOF
}

for row in "${SYSTEMS[@]}"; do
  read -r sys nuc sshort slong docpu <<< "$row"
  src="$PERF/$sys"
  conf="$src/init_conf_${sys}.dat"
  [ -f "$conf" ] || { echo "# missing $conf - skipping $sys" >&2; continue; }

  for cfg in \
     "cpu|backend = CPU|$OXC" \
     "float_verlet|backend = Metal\nMetal_avoid_cpu_calculations = 1\nbackend_precision = float\nMetal_list = verlet|$OXM" \
     "float_edge|backend = Metal\nMetal_avoid_cpu_calculations = 1\nbackend_precision = float\nMetal_list = edge|$OXM" \
     "mixed_verlet|backend = Metal\nMetal_avoid_cpu_calculations = 1\nbackend_precision = mixed\nMetal_list = verlet|$OXM" \
     "mixed_edge|backend = Metal\nMetal_avoid_cpu_calculations = 1\nbackend_precision = mixed\nMetal_list = edge|$OXM" \
     "hardmixed|backend = Metal\nMetal_avoid_cpu_calculations = 1\nbackend_precision = hardmixed\nMetal_list = verlet|$OXO" \
  ; do
     IFS="|" read -r name bk bin <<< "$cfg"
     { [ "$name" = "cpu" ] && { [ "$docpu" = "0" ] || [ -z "$bin" ]; }; } && continue
     d="$WORK/${sys}_${name}"; mkdir -p "$d"
     cp "$conf" "$d/${sys}_conf.dat"; cp "$src"/topology_${sys}.top "$d/"
     pfx=""; [ "$name" = "hardmixed" ] && pfx="OMP_NUM_THREADS=${OMP_NUM_THREADS:-6}"

     bs=99999; bl=99999; U=""
     for rep in 1 2; do
        mkinput "$d" "$(printf '%b' "$bk")" "$sshort" "$sys"
        ws=$( cd "$d" && env $pfx /usr/bin/time -p "$bin" input 2>&1 >/dev/null | awk '/^real/{print $2}' )
        mkinput "$d" "$(printf '%b' "$bk")" "$slong" "$sys"
        wl=$( cd "$d" && env $pfx /usr/bin/time -p "$bin" input 2>&1 >/dev/null | awk '/^real/{print $2}' )
        U=$(tail -1 "$d/e.dat" 2>/dev/null | awk '{print $4}')
        bs=$(awk -v a="$ws" -v b="$bs" 'BEGIN{print (a<b)?a:b}')
        bl=$(awk -v a="$wl" -v b="$bl" 'BEGIN{print (a<b)?a:b}')
     done
     msps=$(awk -v s="$bs" -v l="$bl" -v ns="$sshort" -v nl="$slong" 'BEGIN{printf "%.4f",(l-s)/(nl-ns)*1000.0}')
     tot=$(awk -v m="$msps" 'BEGIN{printf "%.3f", m*1e6/1000.0/3600.0}')
     echo "$sys,$nuc,$name,$sshort,$slong,$bs,$bl,$msps,$tot,$U"
  done
done
