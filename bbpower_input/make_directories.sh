#!/bin/bash

name_run=230324/
base_folder=/global/cfs/cdirs/act/data/iabril/BBPower/

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