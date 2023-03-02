#!/bin/bash

base_folder=/global/cfs/cdirs/act/data/iabril/BBPower/230227_d10
#input_folder=/global/common/software/act/iabril/python/DustNGwBBP/bbpower_input

machine=perl
experiment=so
nside=256
cov_corr=fsky
lmin=30
nbands=9
dell=30
dellnmt=10
pol=B
weight=Cl

mtype=full
apodeg=nan
smooth=20.0
templ=d10
# w_type=w #wt
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
for w_type in $wtypes
do

    name_run=${experiment}_${nside}_${cov_corr}_${templ}_${lmin}_${nbands}_${dell}_${mtype}_${apodeg}_${smooth}_${dellnmt}_${pol}_${ctype}_${weight}
    name_c=${name_run}_${cross}_${moments}_${w_type}_${lminbb}_${lmaxbb}_${bands}
    name_cls=${name_run}_${w_type}
    echo ${base_folder}
    echo ${name_c}

#    if [ -f ${base_folder}/chains/${name_c}.npz.h5 ]; then
 #       rm ${base_folder}/chains/${name_c}.npz*
 #   fi

    # Run pipeline
    python -m bbpower BBCompSep --cells_coadded=$base_folder/${name_cls}${all}  --cells_noise=$base_folder/${name_cls}${noise} --cells_fiducial=$base_folder/${name_cls}${fid}    --config=$base_folder/config_files/${name_c}.yml      --param_chains=$base_folder/chains/${name_c}.npz       --config_copy=$base_folder/temp/config_copy.yml

done
#done
#done
#done
#done