1) cgna_SEQave_gs_20bp_17jx_noends.npy and cgna_SEQave_cov_20bp_17jx_noends.npy contain the ground state and full covariance matrix for 20bp sequences.
gs and cov are averaged over 10000 sequences. only 17 junctions are considered (i.e ends are excluded) and the phosphates coordinates have been integrated out

2) cgna_ave_gs.npy and cgna_ave_cov.npy are the average mean and covariance for a given junction predicted by cgna. They were obtained from the cgna_SEQave_gs_20bp_17jx_noends.npy and cgna_SEQave_cov_20bp_17jx_noends.npy by averaging over the central (farther from the ends) 6 bps.

3) to load these files, use np.load(file)
