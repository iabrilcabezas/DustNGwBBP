## --cell-noise and fiducial can be rubbish because i am not using them

base_folder=/global/cfs/cdirs/act/data/iabril/BBPower/bestfit_ptes/
rubbish=so_256_w2_p353_30_9_30_soflat_5.0_0.4_10_B_dcs_Cl_w_fid.fits
config_file=/global/common/software/act/iabril/python/DustNGwBBP/pte/config_bestfitpte.yml
name_run=so_256_w2_p353_30_9_30_soflat_5.0_0.4_10_B_dcs_Cl
wtypes='w dfwt'

for no_sim in {999..4000}
do
echo ${no_sim}
for w_type in $wtypes
do
python -m bbpower BBCompSep --cells_fiducial=${base_folder}/${rubbish} --cells_noise=${base_folder}/${rubbish} --cells_coadded=${base_folder}/${no_sim}/${w_type}/${name_run}_${w_type}_tot.fits --config=${config_file} --output_dir=${base_folder}/${no_sim}/${w_type}/
done
done