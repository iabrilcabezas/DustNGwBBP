# get mask:
w = get_mask_full(NSIDE)# get_mask_bicep(NSIDE,  apo_deg )
# get bandpowers
b, ell_eff = get_bandpowers(NSIDE, width_l)
# get template:
t = get_template(machine, NSIDE, apo_deg , smooth_deg)

f_d2 = nmt.NmtField(w * w, None, spin = 0)
f_dtilde2 = nmt.NmtField((w * t)**2, None, spin = 0)

mw2_matrix =  get_couplingmatrix(f_d2, f_d2, b)
mwt2_matrix = get_couplingmatrix(f_dtilde2, f_dtilde2, b)

np.savetxt(output_path + '_'.join([mtype, name_run]) + '_couplingM_w.txt', mw2_matrix)
np.savetxt(output_path + '_'.join([mtype, name_run]) + '_couplingM_wt.txt', mwt2_matrix)