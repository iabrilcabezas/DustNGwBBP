#!/bin/bash

base_folder=/global/cfs/cdirs/act/data/iabril/BBPower/pipeline_test
#input_folder=/global/common/software/act/iabril/python/DustNGwBBP/bbpower_input

w_type=w
experiment=so
nside=256
lmin=2
nbands=50
dell=10
smooth=40.0

test_type=${experiment}_full/$w_type
config_type=bbpw_${experiment}_full_${w_type}
name_run=perl_${experiment}_${nside}_${lmin}_${nbands}_${dell}_full_nan_${smooth}_10_B

name_config=${name_run}_all_Cl_${w_type}
name_chain=${experiment}_${nside}_${lmin}_${nbands}_${dell}_${smooth}_${w_type}

all=_all_Cl_${w_type}_tot.fits
noise=_all_Cl_${w_type}_noi.fits 
fid=_all_Cl_${w_type}_fid.fits

# rm $base_folder/$test_type/param_chains*

# Run pipeline
python -m bbpower BBCompSep --cells_coadded=$base_folder/$name_run$all  --cells_noise=$base_folder/$name_run$noise --cells_fiducial=$base_folder/$name_run$fid    --config=$base_folder/config_files/${name_config}.yml      --param_chains=$base_folder/chains/${name_chain}.npz       --config_copy=$base_folder/temp/config_copy.yml

# This plots the results of the pipeline
# i'll do the plotting myself, thank you
# python -m bbpower BBPlotter  --cells_coadded_total=$base_folder/$name_run$all  --cells_coadded=$base_folder/$name_run$all   --cells_null=$base_folder/$name_run$all     --cells_noise=$base_folder/$name_run$noi        --cells_fiducial=$base_folder/$name_run$fid           --param_chains=$base_folder/chains/${name_chain}.npz            --plots=$base_folder/temp/plots                  --plots_page=$base_folder/temp/plots/plots_page.html                     --config=$base_folder/config_files/${name_config}.yml

# # Check the final plots exist
# if [ ! -f $base_folder/temp/plots/triangle.png ]; then
#     echo "Test did not pass"
# else
#     echo "Test passed"
# fi