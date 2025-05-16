#!/bin/bash

#SBATCH --array=1-800
#SBATCH -N 1 # node count
#SBATCH -c 1
#SBATCH -t 48:00:0
#SBATCH --mem=10000
#SBATCH --mail-type=ALL
#SBATCH -o logs/no_osm_sweep_Keq_1_no_sat/%a.out

#Example script running simulation sweep of linear model with Keq=1
table_file='../parameters/hpc_tables/no_osm_sweep_Keq_1_no_sat.csv'
results_folder='results/no_osm_sweep_Keq_1_no_sat'
mkdir $results_folder

true_id=$((SLURM_ARRAY_TASK_ID + 0))

ml matlab

matlab -nojvm -r "addpath(genpath('waste_crossfeeding_evolution')); chemostat_invasion_batch_run( '$table_file' , $true_id , '$results_folder',50); quit();"

