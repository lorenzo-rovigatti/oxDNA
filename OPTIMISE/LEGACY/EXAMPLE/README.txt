To run the example:

1) copy config.txt and the two other files to a directory.
2) enter that directory and run 

	bash path/to/opti/optimise.sh opti_input_example.txt > OptiLog
	
3) > OptiLog is not necessary, but recommended

The example runs a very short test optimisation (2 cycles, 10 rew steps, 2 replicas, 1e6 timesteps).
of the average propeller, rise, roll and twist (1,7,8,11) by using a few parameters (see last section of the input script).
MODE ave = target SD AVERAGED (see README in Params folder) averages (and covs).
IDS_COV is commented out because we are not targeting any covariance entry.  
