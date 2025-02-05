#!/bin/sh
#
#SBATCH --job-name="run_test"
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --account=innovation

source ./modules.sh
/home/fgreco/progs/julia-1.11.3/bin/julia --project=../ -O0 --check-bounds=no -e 'include("../scripts/main.jl")'
