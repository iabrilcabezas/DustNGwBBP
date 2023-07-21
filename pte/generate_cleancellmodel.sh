# module load python
# conda activate /global/cfs/cdirs/act/software/iabril/condaenvs/BBenv
# cd /global/cfs/cdirs/act/software/iabril/condaenvs/github_reps/BBPower

# ### FIRST PART, GENERATE DATA
# # --cells_* can be rubbish

base_folder=/global/cfs/cdirs/act/data/iabril/BBPower/230414
rubbish=so_256_w2_p353_30_9_30_soflat_5.0_0.4_10_B_dcs_Cl_w


out_base=/global/cfs/cdirs/act/data/iabril/BBPower/230525_sample/sims/ #03_pte

wtypes='w dfwt'

for i in {10000..29999}
do
echo ${i}
for w_type in $wtypes
do

outdir=${out_base}/${i}/${w_type}

python -m bbpower BBCompSep --cells_fiducial=${base_folder}/${rubbish}_fid.fits --cells_noise=${base_folder}/${rubbish}_noi.fits --cells_coadded=${base_folder}/${rubbish}_tot.fits --cells_coadded_cov=${base_folder}/${rubbish}_tot.fits --config=${outdir}/bestfitcell.yml --output_dir=${outdir}

done
done
