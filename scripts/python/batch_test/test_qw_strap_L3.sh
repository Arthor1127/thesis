#!/bin/bash
#$ -cwd
#$ -N qw_strap_2L3
#$ -j y
#$ -S /bin/bash
#$ -q copahue,highmem
#$ -l memoria_a_usar=24G
#$ -pe neworte 32
#$ -V

# Single eigenmode run: quarter-wave resonator, strap at x0 = 2L/3.
#
# x0 measured from the SHORT end (sg1), matching how
# RouteMeanderGrounded's 'positions' and qw_sweep_plan.py both measure
# it. x0 = 2*5.92905/3 = 3.9527mm.
#
# 2L/3 is the interesting point: BOTH segments resonate at 3*f1
# simultaneously (a degenerate pair), because it is the 3*lambda/4
# mode's own node, so the strap barely perturbs that mode.
#   segment A [0, x0]  short-short -> v/(2*x0)     = 15.000 GHz
#   segment B [x0, L]  short-open  -> v/(4*(L-x0)) = 15.000 GHz
# --starting-freq/--target-freq-ghz center on 15GHz. NOTE: an earlier
# version of this script targeted 7.5GHz, the prediction for L/3 (strap
# measured from the OPEN end) - wrong for x0 measured from the short
# end, and ~7.5GHz away from the real mode. That is well clear
# of the box mode confirmed at ~12.3-12.4GHz (chip-resize diagnostic:
# stayed put in frequency-normalized terms while 3*f1 moved as expected
# for a real resonator mode - see boxdiag results), so there is little
# risk --target-freq-ghz accidentally selects the box mode here.
# --number-of-freqs 6 still returns enough neighbours to see the box
# mode explicitly in the output as a sanity check.
#
# NODE CHOICE: peak memory for this design has run ~10G in every prior
# eigenmode test, comfortably under the 24G ceiling requested - no need
# for highmem/compute-4-18's 503G, and after the RAM-courtesy point
# raised this session, using copahue's ordinary nodes for a small job is
# the better citizen choice, leaving compute-4-18 free for jobs that
# actually need it.
#
# "As many cores as possible" was read as "most PHYSICAL cores", not
# highest NCPU: compute-4-3/4/5 show NCPU=40 but NCOR=20 (2x
# hyperthreading) - oversubscribing MPI ranks past physical core count
# tends to hurt on this workload, per the same overhead reasoning behind
# the rank-scaling experiment, so those are not actually "more cores"
# for this purpose. Among genuine 32-physical-core nodes: compute-4-17
# was at load 89.5 (something else is hammering it) and compute-4-15 at
# load 8 - both excluded. compute-4-14 and compute-4-16 were idle
# (load 0.01) at submission time - restricted to exactly these two by
# EXCLUDING every other host, the '!host&' chain syntax already
# confirmed working in stageB_driven.sh this session (an OR-based
# inclusion syntax was tried first and reverted here - not worth
# risking on an operator this project has not actually confirmed works).

conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH

export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg

PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/batch_test"
cd "$SCRIPT_DIR"

OUTDIR=$HOME/sweep_qw/geom_L3strap
mkdir -p $OUTDIR

echo "=== Quarter-wave, strap at x0=2L/3=3.9527mm (from SHORT end) ==="
echo "=== Predicted mode: 15.000 GHz (degenerate pair) (box mode confirmed separately at ~12.3GHz) ==="
echo "=== Running on host: $(hostname) with $NSLOTS ranks ==="

python build_and_eigenmode_qw.py \
    --ground-pos 3.9527 \
    --outdir $OUTDIR \
    --palace-bin "$PALACE_BIN" \
    --num-cpus $NSLOTS \
    --total-length-mm 5.92905 \
    --starting-freq 15.0e9 \
    --number-of-freqs 6 \
    --bg-min-size-um 8 --bg-max-size-um 100 \
    --target-freq-ghz 15.0
