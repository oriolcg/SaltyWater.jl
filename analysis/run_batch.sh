#!/bin/sh
#
#SBATCH --job-name="run_test"
#SBATCH --partition=compute
#SBATCH --time=00:01:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --account=innovation

source ./modules.sh

Pressures=("5000" "6000" "7000" "8000")
Velocities=("0.09" "0.129" "0.258")
for p in ${Pressures[@]}
do
    for v in ${Velocities[@]}
    do
        echo $p $v
    done
done

#julia --project=../ -O0 --check-bounds=no -e 'include("../scripts/main.jl")'
