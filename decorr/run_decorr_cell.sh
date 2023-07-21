#(BBenv)iabril@nersc:login22 [04:28:57] [/global/cfs/cdirs/act/software/iabril/condaenvs/github_reps/BBPower]

# module load python
# conda activate /global/cfs/cdirs/act/software/iabril/condaenvs/BBenv
# cd /global/cfs/cdirs/act/software/iabril/condaenvs/github_reps/BBPower

# ### FIRST PART, GENERATE DATA WITH DECORR
# # --cells_* can be rubbish
# cd /global/cfs/cdirs/act/software/iabril/condaenvs/github_reps/BBPower

# base_folder=/global/cfs/cdirs/act/data/iabril/BBPower/230414
# rubbish=so_256_w2_p353_30_9_30_soflat_5.0_0.4_10_B_dcs_Cl_w

# config_file=/global/common/software/act/iabril/python/DustNGwBBP/decorr/bbcompsep_decorr.yml
# outdir=/global/cfs/cdirs/act/data/iabril/BBPower/230524_decorr
# python -m bbpower BBCompSep --cells_fiducial=${base_folder}/${rubbish}_fid.fits --cells_noise=${base_folder}/${rubbish}_noi.fits --cells_coadded=${base_folder}/${rubbish}_tot.fits --cells_coadded_cov=${base_folder}/${rubbish}_tot.fits --config=${config_file} --config_copy=${outdir}/temp/config_copy_generateD.yml --output_dir=${outdir} 

# add covariance to cells:
#cd /global/common/software/act/iabril/python/DustNGwBBP/
#python run_decorr.py
## SECOND PART, RUN BBCOMPSEP:

base_folder=/global/cfs/cdirs/act/data/iabril/BBPower/230414
base_decorr=/global/cfs/cdirs/act/data/iabril/BBPower/230622_decorr

all=_tot.fits
noise=_noi.fits 
fid=_fid.fits

# wtypes='w dfwt' # ' wt '
w_type='dfwt'

name_run=so_256_w2_p353_30_9_30_soflat_5.0_0.4_10_B_dcs_Cl
config_file=/global/common/software/act/iabril/python/DustNGwBBP/decorr/bbcompsep_decorr.yml

#for w_type in $wtypes
#do

   name_cls=${name_run}_${w_type}

     if [ ! -d ${base_decorr}/chains/${w_type} ]; then
             mkdir ${base_decorr}/chains/${w_type}
     fi

    python -m bbpower BBCompSep --cells_coadded=$base_decorr/${name_cls}${all}  --cells_coadded_cov=$base_decorr/${name_cls}${all}  --config=${config_file}      --output_dir=$base_decorr/chains/${w_type}/    --config_copy=$base_decorr/temp/config_copy_${w_type}.yml --cells_noise=$base_folder/${name_cls}${noise} --cells_fiducial=$base_decorr/cells_model.fits # (no needed for chi2)

#done