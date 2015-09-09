### Hooks for the editor to set the default target
current: target

target pngtarget pdftarget vtarget acrtarget: NIH_example.Rout 

##################################################################

# make files

Sources = Makefile .gitignore 

##################################################################

# compiler.mk is intended to be a local file; it uses clang by default, but you can just edit and save it on any particular machine. Don't add it to sources

include compiler.mk
compiler.mk:
	/bin/cp clang.mk $@

# Add your own compiler makefile to this list to share it between machines
Sources += clang.mk gcc.mk

# Old C stuff that we may not need
Sources += exe.mk

######################################################################

# Make output file from fake parameters

format_functions = run_sherif_FCT.R format_NIH_FCT.R format_NIH.R
Sources += $(format_functions)
format_functions.Rout: $(format_functions)
	$(run-R)

Sources += run_sherif_FCT.R Calibration/loadParam_FCT.R

Sources += NIH_example.R
Sources += param_model.csv param_simul.csv prediction_date.csv

NIH_example.Rout: format_functions.Rout Calibration/loadParam_FCT.Rout
NIH_example.Rout: param_model.csv param_simul.csv prediction_date.csv
NIH_example.Rout: NIH_example.R
	-/bin/rm -rf NIH_example_dir.old
	-/bin/mv -f NIH_example_dir NIH_example_dir.old
	mkdir NIH_example_dir
	$(run-R)

NIH_example_dir: NIH_example.Rout ; 


### Temporary copying rules ###

NIH_example.R: ~/Dropbox/SHERIF/run_sherif.R
	/bin/cp $< $@
%.R: ~/Dropbox/SHERIF/%.R
	/bin/cp $< $@

%.csv: ~/Dropbox/SHERIF/%.csv
	/bin/cp $< $@

Calibration/%.R: ~/Dropbox/SHERIF/Calibration/%.R
	/bin/cp $< $@

# JD git rules

ms = ../makestuff
-include $(ms)/git.mk
-include $(ms)/RR.mk
