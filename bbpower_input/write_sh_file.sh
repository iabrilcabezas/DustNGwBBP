#!/bin/bash

base_folder=/global/cfs/cdirs/act/data/iabril/BBPower/pipeline_test
#input_folder=/global/common/software/act/iabril/python/DustNGwBBP/bbpower_input

machine=cori
experiment=so
nside=256
lmin=2
nbands=50
dell=10
mtype=full
apodeg=nan
smooth=20.0
dellnmt=10
pol=B
weight=Cl

ctype=dcs
#w_type=wt
cross=C
moments=M

all=_tot.fits
noise=_noi.fits 
fid=_fid.fits

name_run=${machine}_${experiment}_${nside}_${lmin}_${nbands}_${dell}_${mtype}_${apodeg}_${smooth}_${dellnmt}_${pol}_${ctype}_${weight}

wtypes='w wt'

for w_type in $wtypes
do
    name_c=${name_run}_${cross}_${moments}_${w_type}
    name_cls=${name_run}_${w_type}

    if [ -f $base_folder/chains/${name_c}.npz* ]; then
        rm $base_folder/chains/${name_c}.npz*
    fi

    # Run pipeline
    python -m bbpower BBCompSep --cells_coadded=$base_folder/${name_cls}${all}  --cells_noise=$base_folder/${name_cls}${noise} --cells_fiducial=$base_folder/${name_cls}${fid}    --config=$base_folder/config_files/${name_c}.yml      --param_chains=$base_folder/chains/${name_c}.npz       --config_copy=$base_folder/temp/config_copy.yml
done
# This plots the results of the pipeline
# i'll do the plotting myself, thank you
# python -m bbpower BBPlotter  --cells_coadded_total=$base_folder/$name_run$all  --cells_coadded=$base_folder/$name_run$all   --cells_null=$base_folder/$name_run$all     --cells_noise=$base_folder/$name_run$noi        --cells_fiducial=$base_folder/$name_run$fid           --param_chains=$base_folder/chains/${name_chain}.npz            --plots=$base_folder/temp/plots                  --plots_page=$base_folder/temp/plots/plots_page.html                     --config=$base_folder/config_files/${name_config}.yml

# # Check the final plots exist
# if [ ! -f $base_folder/temp/plots/triangle.png ]; then
#     echo "Test did not pass"
# else
#     echo "Test passed"
# fi