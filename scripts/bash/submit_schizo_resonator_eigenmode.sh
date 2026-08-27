#!/bin/bash
# ======================================================
#  Cluster launcher for schizo_resonator_eigenmode.py (SQDMetal + Palace)
#  Compatible with Copahue/Caulle/CPUinGPU queues
# ======================================================
# Palace parallelizes internally via MPI (it is invoked as
# `palace -np NUM_CPUS -nt 1 config.json` by SQDMetal), so this
# is a single SGE job that reserves NUM_CPUS slots on one node
# and lets Palace use all of them - no job array needed here.
#
# IMPORTANT: fill in PALACE_DIR below with the real path from:
#   spack location -i palace
# once your `spack install palace +superlu-dist +slepc %gcc@16.1.0` finishes.
#
# Usage:
#   ./submit_schizo_resonator_eigenmode.sh <queue_name> <num_cpus> <output_folder>
# Example:
#   ./submit_schizo_resonator_eigenmode.sh copahue 32 schizo_eigen_2026_08_12
# ======================================================

if [ $# -ne 3 ]; then
  echo "Usage: $0 <queue_name> <num_cpus> <output_folder>"
  exit 1
fi

QUEUE_NAME=$1
NUM_CPUS=$2
OUTPUT=$3

HOME_DIR=/home/arturo.ruiz/quantum_circuits/schizo_resonators
PYFILE=schizo_resonator_eigenmode.py
PALACE_DIR=/home/arturo.ruiz/spack/opt/spack/linux-x86_64/palace-XXXX/bin/palace   # <-- fill in from `spack location -i palace`
OUTPUT_DIR=${HOME_DIR}/${OUTPUT}
mkdir -p "$OUTPUT_DIR"

echo "Submitting schizo_resonator_eigenmode Palace job"
echo "Queue: $QUEUE_NAME"
echo "Cores requested (matches Palace -np): $NUM_CPUS"
echo "Output base directory: ${OUTPUT_DIR}"
echo

# ======================================================
#  Create master execution script ("script-var")
# ======================================================
cat > script-var <<EOF
#!/bin/bash
#$ -clear
#$ -S /bin/bash
#$ -N schizo_eigen
#$ -cwd
#$ -pe mpi ${NUM_CPUS}
#$ -l mem=4G
#$ -j y
#$ -o output_log
#$ -q ${QUEUE_NAME}

export OMP_NUM_THREADS=1
export QT_QPA_PLATFORM=offscreen

# --- Environment setup ---
TMPDIR=/tmp/\$USER/\$JOB_ID
mkdir -p \$TMPDIR

source /share/apps/miniconda3/etc/profile.d/conda.sh
conda activate quantum

# --- Directories ---
HOME_DIR=${HOME_DIR}
OUTPUT_DIR=${OUTPUT_DIR}
JOB_DIR=\$OUTPUT_DIR/job_\${JOB_ID}
mkdir -p \$JOB_DIR

# Copy Python code for reproducibility
cp \$HOME_DIR/${PYFILE} \$JOB_DIR/
cd \$JOB_DIR || exit 1

# --- Log job info ---
{
  echo "Running on \$(hostname)"
  echo "Slots (NSLOTS): \$NSLOTS"
  echo "Python: \$(which python)"
  echo "Conda env: \$(conda info --envs | grep '*')"
} > job_log.txt

# --- Run simulation ---
python ${PYFILE} --num-cpus \$NSLOTS --palace-dir ${PALACE_DIR} \\
    --sim-dir res_eigen_test --outdir \$JOB_DIR \\
    --total-length 3.6mm --ground-strap-positions 1.2mm 2.4mm --ground-strap-width 10um \\
    --starting-freq 10e9 --number-of-freqs 10 --solver-tol 1.0e-7 \\
    > output_log.txt 2>&1

# --- Copy results to final directory ---
cp -r \$JOB_DIR \$OUTPUT_DIR/

# --- Clean temporary directory ---
rm -rf \$TMPDIR
EOF

# ======================================================
#  Copy launcher components to output directory
# ======================================================
cp -p "$PYFILE" script-var "$OUTPUT_DIR" 2>/dev/null
cd "$OUTPUT_DIR" || exit 1

# ======================================================
#  Submit job to cluster
# ======================================================
qsub script-var

echo
echo "Submitted 1 job requesting ${NUM_CPUS} cores on queue ${QUEUE_NAME}."
