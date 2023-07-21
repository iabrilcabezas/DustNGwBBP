#!/bin/bash

base_folder=/global/cfs/cdirs/act/data/iabril/BBPower/230414
#input_folder=/global/common/software/act/iabril/python/DustNGwBBP/bbpower_input

machine=perl
experiment=so
nside=256
cov_corr=w2
lmin=30
nbands=9
dell=30
dellnmt=10
pol=B
weight=Cl

mtype=soflat
apodeg=5.0
smooth=0.4
templ=p353
# w_type=dfwt # 00 #w #wt
ctype=dcs
cross=C
moments=M
decor=0

bands=all #MF1_UHF2 # LF1_LF2_MF1_MF2_UHF2 
lminbb=30
lmaxbb=300

all=_tot.fits
noise=_noi.fits 
fid=_fid.fits

wtypes='w dfwt' # ' wt '
ctypes='dc0 dcs'
mmts='0 M'
crosses='0 C'
decorrs='0 D'
smoothes=( "0.4" "10.0")

#for ctype in $ctypes
#do
#for smooth in "${smoothes[@]}"
#do
# for cross in $crosses
# do
#for moments in $mmts
#do
#for decor in $decorrs
#do
for w_type in $wtypes
do

    name_run=${experiment}_${nside}_${cov_corr}_${templ}_${lmin}_${nbands}_${dell}_${mtype}_${apodeg}_${smooth}_${dellnmt}_${pol}_${ctype}_${weight}
    name_c=${name_run}_${cross}_${moments}_${decor}_${w_type}_${lminbb}_${lmaxbb}_${bands}
    name_cls=${name_run}_${w_type}
    echo ${base_folder}
    echo ${name_c}

   # if [ -f ${base_folder}/chains/${name_c}.npz.h5 ]; then
   #         rm ${base_folder}/chains/${name_c}.npz*
   #fi
# new horrible bbpower:
    if [ ! -d ${base_folder}/chains/${name_c} ]; then
            mkdir ${base_folder}/chains/${name_c}
    fi

    # Run pipeline
    python -m bbpower BBCompSep --cells_coadded=$base_folder/${name_cls}${all} --cells_coadded_cov=$base_folder/${name_cls}${all} --cells_noise=$base_folder/${name_cls}${noise} --cells_fiducial=$base_folder/${name_cls}${fid}    --config=$base_folder/config_files/${name_c}.yml  --output_dir=${base_folder}/chains/${name_c}    --config_copy=$base_folder/temp/config_copy.yml #   --param_chains=$base_folder/chains/${name_c}.npz   
done
# done
#done
#done
#done
#done
#done