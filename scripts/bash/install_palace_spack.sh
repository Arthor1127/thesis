#!/bin/bash
# ======================================================
#  Build Palace via Spack as an SGE job (not on frontend)
# ======================================================
# Prereqs (run once, on the frontend, needs internet):
#   module load miniconda
#   git clone -c feature.manyFiles=true https://github.com/spack/spack.git ~/spack
#   . ~/spack/share/spack/setup-env.sh
#   spack compiler find
#   spack external find openmpi   # or mpich - reuse cluster's MPI/libs
#   spack fetch palace            # downloads sources while on the frontend
#
# Then submit this script:
#   qsub install_palace_spack.sh
# ======================================================
#$ -clear
#$ -S /bin/bash
#$ -N spack_install_palace
#$ -cwd
#$ -q copahue
#$ -pe mpi 8
#$ -l mem=4G
#$ -l h_rt=12:00:00
#$ -j y
#$ -o install_palace.log

. ~/spack/share/spack/setup-env.sh

echo "Building on $(hostname) with $NSLOTS build jobs"
spack install -j "$NSLOTS" palace

echo "Install location:"
spack location -i palace
