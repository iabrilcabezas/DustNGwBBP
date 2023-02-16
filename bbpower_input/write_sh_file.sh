#!/bin/bash

base_folder=/global/cfs/cdirs/act/data/iabril/BBPower/230215_now2factor
#input_folder=/global/common/software/act/iabril/python/DustNGwBBP/bbpower_input

machine=perl
experiment=so
nside=256
cov_corr=None
lmin=2
nbands=50
dell=10
dellnmt=10
pol=B
weight=Cl

mtype=so
apodeg=5.0
smooth=20.0
#w_type=wt
ctype=dcs
cross=0
moments=0

bands=all #MF1_UHF2 #
lminbb=30
lmaxbb=300

all=_tot.fits
noise=_noi.fits 
fid=_fid.fits

wtypes='w wt'
ctypes='dc0 dcs'
mmts='0 M'
crosses='0 C'
smoothes=( "10.0" "20.0" "40.0")

#for ctype in $ctypes
#do
#for smooth in "${smoothes[@]}"
#do
# for cross in $crosses
# do
#for moments in $mmts
#do
#for w_type in $wtypes
#do

w_type=w

name_run=${experiment}_${nside}_${cov_corr}_${lmin}_${nbands}_${dell}_${mtype}_${apodeg}_${smooth}_${dellnmt}_${pol}_${ctype}_${weight}
name_c=${name_run}_${cross}_${moments}_${w_type}_${lminbb}_${lmaxbb}_${bands}
name_cls=${name_run}_${w_type}
echo ${name_c}

if [ -f ${base_folder}/chains/${name_c}.npz.h5 ]; then
    rm ${base_folder}/chains/${name_c}.npz*
fi

# Run pipeline
python -m bbpower BBCompSep --cells_coadded=$base_folder/${name_cls}${all}  --cells_noise=$base_folder/${name_cls}${noise} --cells_fiducial=$base_folder/${name_cls}${fid}    --config=$base_folder/config_files/${name_c}.yml      --param_chains=$base_folder/chains/${name_c}.npz       --config_copy=$base_folder/temp/config_copy.yml
#done
#done
#done
#done
#done
# This plots the results of the pipeline
# i'll do the plotting myself, thank you
# python -m bbpower BBPlotter  --cells_coadded_total=$base_folder/$name_run$all  --cells_coadded=$base_folder/$name_run$all   --cells_null=$base_folder/$name_run$all     --cells_noise=$base_folder/$name_run$noi        --cells_fiducial=$base_folder/$name_run$fid           --param_chains=$base_folder/chains/${name_chain}.npz            --plots=$base_folder/temp/plots                  --plots_page=$base_folder/temp/plots/plots_page.html                     --config=$base_folder/config_files/${name_config}.yml

# # Check the final plots exist
# if [ ! -f $base_folder/temp/plots/triangle.png ]; then
#     echo "Test did not pass"
# else
#     echo "Test passed"
# fi