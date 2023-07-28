## run bestfit BBCompSep to generate chi2.npz

# --cells_* can be rubbish
#python -m bbpower BBCompSep --cells_fiducial=/global/cfs/cdirs/act/data/iabril/BBPower/230414/so_256_w2_p353_30_9_30_soflat_5.0_0.4_10_B_dcs_Cl_w_fid.fits --cells_noise=/global/cfs/cdirs/act/data/iabril/BBPower/230414/so_256_w2_p353_30_9_30_soflat_5.0_0.4_10_B_dcs_Cl_w_noi.fits --cells_coadded=/global/cfs/cdirs/act/data/iabril/BBPower/230414/so_256_w2_p353_30_9_30_soflat_5.0_0.4_10_B_dcs_Cl_w_tot.fits  --config=/global/common/software/act/iabril/python/DustNGwBBP/pte/config_bestfitpte.yml --output_dir=/global/cfs/cdirs/act/data/iabril/BBPower/230419_PTE/

## --cell-noise and fiducial can be rubbish because i am not using them

base_folder=/global/cfs/cdirs/act/data/iabril/BBPower/230725/sims/
rubbish=so_256_w2_p353_30_9_30_soflat_5.0_0.4_10_B_dcs_Cl_w_fid.fits
config_file=/global/common/software/act/iabril/python/DustNGwBBP/pte/config_bestfitpte
name_run=so_256_w2_p353_30_9_30_soflat_5.0_0.4_10_B_dcs_Cl
bands='LF MF UHF'
# wtypes='w dfwt'
w_type=dfwt

for no_sim in {50000..54999}
do

echo ${no_sim}

python -m bbpower BBCompSep --cells_fiducial=${base_folder}/${rubbish} --cells_noise=${base_folder}/${rubbish} --cells_coadded=${base_folder}/${no_sim}/w/${name_run}_w_tot.fits --config=${config_file}_all.yml --output_dir=${base_folder}/${no_sim}/w/
python -m bbpower BBCompSep --cells_fiducial=${base_folder}/${rubbish} --cells_noise=${base_folder}/${rubbish} --cells_coadded=${base_folder}/${no_sim}/${w_type}/${name_run}_${w_type}_tot.fits --config=${config_file}_all.yml --output_dir=${base_folder}/${no_sim}/${w_type}/

# for band in $bands
# do
# python -m bbpower BBCompSep --cells_fiducial=${base_folder}/${rubbish} --cells_noise=${base_folder}/${rubbish} --cells_coadded=${base_folder}/${no_sim}/w/${name_run}_w_tot.fits --cells_coadded_cov=${base_folder}/${no_sim}/w/${name_run}_w_tot.fits --config=${config_file}_${band}.yml --output_dir=${base_folder}/${no_sim}/w/${band}/
# done
done