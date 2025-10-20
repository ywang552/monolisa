#!/bin/bash
#SBATCH -J MonoLisa
#SBATCH -o /home/ywang552/monolisa/slurm/out/%x-%j.out
#SBATCH -e /home/ywang552/monolisa/slurm/err/%x-%j.err
#SBATCH --cpus-per-task=1
#SBATCH --partition=longjobs

set -euo pipefail

JULIA_EXEC=/home/ywang552/julia-1.9.1/bin/julia
export OMP_NUM_THREADS=1
export JULIA_NUM_THREADS=1

echo ">>> MonoLisa start $(date)"

"$JULIA_EXEC" --project=. main.jl

echo "<<< MonoLisa done  $(date)"
