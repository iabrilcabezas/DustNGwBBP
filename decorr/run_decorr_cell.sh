#(BBenv)iabril@nersc:login22 [04:28:57] [/global/cfs/cdirs/act/software/iabril/condaenvs/github_reps/BBPower]
# --cells_* can be rubbish
#python -m bbpower BBCompSep --cells_fiducial=/global/cfs/cdirs/act/data/iabril/BBPower/230414/so_256_w2_p353_30_9_30_soflat_5.0_0.4_10_B_dcs_Cl_w_fid.fits --cells_noise=/global/cfs/cdirs/act/data/iabril/BBPower/230414/so_256_w2_p353_30_9_30_soflat_5.0_0.4_10_B_dcs_Cl_w_noi.fits --cells_coadded=/global/cfs/cdirs/act/data/iabril/BBPower/230414/so_256_w2_p353_30_9_30_soflat_5.0_0.4_10_B_dcs_Cl_w_tot.fits  --config=/global/common/software/act/iabril/python/DustNGwBBP/decorr/config_decorr.yml --output_dir=/global/cfs/cdirs/act/data/iabril/BBPower/230417/


## RUN BBCOMPSEP:

base_folder=/global/cfs/cdirs/act/data/iabril/BBPower/230414
base_decorr=/global/cfs/cdirs/act/data/iabril/BBPower/230417

# machine=perl
# experiment=so
# nside=256
# cov_corr=w2
# lmin=30
# nbands=9
# dell=30
# dellnmt=10
# pol=B
# weight=Cl

# mtype=soflat
# apodeg=5.0
# # smooth=0.4
# templ=p353
# # w_type=dfwt # 00 #w #wt
# ctype=dcs
# cross=C
# moments=0
# decor=D

# bands=all
# lminbb=30
# lmaxbb=300

all=_tot.fits
noise=_noi.fits 
fid=_fid.fits

#wtypes='w dfwt' # ' wt '

name_run=so_256_w2_p353_30_9_30_soflat_5.0_0.4_10_B_dcs_Cl

#for w_type in $wtypes
#do
w_type='w'

name_cls=${name_run}_${w_type}

python -m bbpower BBCompSep --cells_coadded=$base_decorr/${name_cls}${all}  --cells_noise=$base_folder/${name_cls}${noise} --cells_fiducial=$base_folder/${name_cls}${fid}    --config=/global/common/software/act/iabril/python/DustNGwBBP/decorr/bbcompsep_decorr.yml      --output_dir=$base_decorr/chains/       --config_copy=$base_decorr/temp/config_copy.yml

#done