### Sherif

### Hooks for the editor to set the default target
current: target

target pngtarget pdftarget vtarget acrtarget: sherif_test.Rout 

##################################################################

# JD's make files toolbox:

Sources = Makefile .gitignore 
ms = ../makestuff
-include $(ms)/git.def


##################################################################
###
###   Build the R library for SHERIF model from C++ code (with Rcpp)
###
##################################################################

lib:
	mkdir $@
	-cp -r ../../lib/sherif $@

# individual.cpp defines the class of individuals
# simulator.cpp runs the model
# mc.cpp is a Monte Carlo wrapper for simulator

libroot = EventNumberContainer dcTools RV globalVar individual mc dcDataFrame simulator dcMatrix spatialSim
libcpp = $(libroot:%=%.cpp) Rwrap_sherif.cpp
libh = $(libroot:%=%.h)
Sources += make_library.R $(libcpp) $(libh) Makevars

lib/sherif: 
	$(MAKE) make_library.Rout

make_library.Rout: $(libh) $(libcpp) lib Makevars make_library.R
	-/bin/rm -rf sherif
	$(run-R)
	cp Makevars sherif/src
	R CMD build sherif
	R CMD check sherif
	R CMD INSTALL -l ./lib sherif
	touch lib/sherif
	touch $@

######################################################################
###
### Test if sherif R library works correctly
###
######################################################################

SpatialCSV = distLocations.csv gravity_cst.csv initLocations.csv nLocations.csv param_model.csv param_simul.csv popLocations.csv prediction_date.csv

sherif_test.Rout: sherif_test.R make_library.Rout
	$(run-R)



##################################################################
# JD git rules
-include $(ms)/git.mk
-include $(ms)/visual.mk
-include $(ms)/RR.mk
-include $(ms)/compare.mk
-include $(ms)/local.mk

Makefile: lib ../makestuff

../makestuff: ../%:
	cd .. && git clone git@github.com:dushoff/$*.git
