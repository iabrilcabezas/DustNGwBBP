module load python
conda activate /global/cfs/cdirs/act/software/iabril/condaenvs/BBenv
cd /global/cfs/cdirs/act/software/iabril/condaenvs/github_reps/BBPower

# ### FIRST PART, GENERATE DATA (without  DECORR)
# # --cells_* can be rubbish
cd /global/cfs/cdirs/act/software/iabril/condaenvs/github_reps/BBPower

base_folder=/global/cfs/cdirs/act/data/iabril/BBPower/230414
rubbish=so_256_w2_p353_30_9_30_soflat_5.0_0.4_10_B_dcs_Cl_w

config_file=/global/common/software/act/iabril/python/DustNGwBBP/pte/bbcompsep_baseline.yml
outdir=/global/cfs/cdirs/act/data/iabril/BBPower/230525_sample
python -m bbpower BBCompSep --cells_fiducial=${base_folder}/${rubbish}_fid.fits --cells_noise=${base_folder}/${rubbish}_noi.fits --cells_coadded=${base_folder}/${rubbish}_tot.fits --cells_coadded_cov=${base_folder}/${rubbish}_tot.fits --config=${config_file} --config_copy=${outdir}/temp/config_copy_generateD.yml --output_dir=${outdir} 
