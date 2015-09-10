################################################################################
###
###   MAKEFILE FOR CALIBRATION
###
###   2015-08-20 by David Champredon
###
################################################################################


### Calibration to generated data
### (use summary statistics #1)
###
calibration.out: calibration_synlik.R inference_FCT.R wrap_FCT.R \
generate-target-data.R synlik_stats.R \
param_true4fit.csv generate-target-data-param.csv \
param_simul.csv lhsParam.csv ncpu.csv
	rm -f calib_results_1.pdf  
	Rscript $< 1 > $@ &   # <-- * WARNING * Uses calibration stats of type #1 by default


### Compare the summary statistics for synthetic likelihood.
### The different summary stats functions are defined in 'synlik_stats.R'
### The results of the comparison are analyzed with 'read_compstats.R'
###
compstats: calibration_synlik.R inference_FCT.R wrap_FCT.R generate-target-data.R synlik_stats.R \
param_true4fit.csv generate-target-data-param.csv param_simul.csv lhsParam.csv ncpu.csv \
						read_compstats.R
	rm -f calib_results_*.pdf
	rm -f check_norm_*.pdf
	rm -f calib_results_*.csv
	Rscript $< 1 > compstats_1.out &
	Rscript $< 2 > compstats_2.out &
	#Rscript $< 3 > compstats_3.out &

# Then run "Rscript read_compstats.R" when compstats_X.out are all done.



### Test the calibration on multiple data set differing by the seed only
###
test_multiseed: test_run_calibMultiSeed.R  inference_FCT.R wrap_FCT.R \
generate-target-data.R synlik_stats.R \
param_true4fit.csv generate-target-data-param.csv \
param_simul.csv lhsParam.csv nseed.csv
	Rscript $< &

### Then "Rscript test_run_calibMultiSeed_analyze.R"


### Clean-up
###
clean:
	rm -f calib_results_*.pdf
	rm -f check_norm_*.pdf
	rm -f calib_results_*.csv
	rm -f *.out
	rm -f generated-target*.pdf
	rm -f compstats*.out
	rm -f compstats*.pdf
	
	
