#!/bin/bash
#SBATCH -J MonoLisa_c4_7kp4
#SBATCH -o /home/ywang552/monolisa/slurm/out/%x-%j.out
#SBATCH -e /home/ywang552/monolisa/slurm/err/%x-%j.err
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=longjobs

set -euo pipefail

JULIA_EXEC=/home/ywang552/julia-1.9.1/bin/julia


echo ">>> MonoLisa start $(date)"

"$JULIA_EXEC" --project=. main.jl

echo "<<< MonoLisa done  $(date)"
