#!/bin/bash

name_run=230417
base_folder=/global/cfs/cdirs/act/data/iabril/BBPower/
base_df=/global/cfs/cdirs/act/data/iabril/DustFilaments/
sims_name=sims_230306/

cd ${base_folder}
mkdir ${name_run}
cd ${name_run}
mkdir chains
mkdir config_files
mkdir temp
mkdir results_posteriors
mkdir results_CL
mkdir results_pte
cd results_pte
mkdir Figures

# cd ${base_df}
# cd ${sims_name}
# mkdir analysis_${name_run}
