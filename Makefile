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

# Build the library

lib:
	mkdir $@

lib/sherif: make_library.Rout ;

libroot = EventNumberContainer dcTools RV globalVar individual mc dcDataFrame simulator dcMatrix
libcpp = $(libroot:%=%.cpp) Rwrap_sherif.cpp
libh = $(libroot:%=%.h)
Sources += make_library.R $(libcpp) $(libh) Makevars

make_library.Rout: $(libh) $(libcpp) make_library.R lib Makevars
	-/bin/rm -rf sherif
	$(run-R)
	cp Makevars sherif/src
	R CMD build sherif
	R CMD check sherif
	R CMD INSTALL -l ./lib sherif
	touch lib/sherif

######################################################################

# NIH_example: Make output file from fake parameters

format_functions = run_sherif_FCT.R loadParam_FCT.R
format_functions += format_NIH_FCT.R format_NIH.R
Sources += $(format_functions)
format_functions.Rout: $(format_functions)
	$(run-R)

Sources += NIH_example.R
Sources += param_model.csv param_simul.csv prediction_date.csv

NIH_example.Rout: lib/sherif
NIH_example.Rout: format_functions.Rout loadParam_FCT.Rout
NIH_example.Rout: param_model.csv param_simul.csv prediction_date.csv
NIH_example.Rout: NIH_example.R
	-/bin/rm -rf NIH_example_dir.old
	-/bin/mv -f NIH_example_dir NIH_example_dir.old
	mkdir NIH_example_dir
	$(run-R)

example_plots.Rout: NIH_example.Rout example_plots.R

NIH_example_dir: NIH_example.Rout ; 

### Temporary copying rules ###

files: $(libh) $(libcpp)
$(libh) $(libcpp):
	/bin/cp ~/Dropbox/SHERIF/$@ .

# JD git rules

ms = ../makestuff
-include $(ms)/git.mk
-include $(ms)/RR.mk

Makefile: ../makestuff

../makestuff: ../%:
	cd .. && git clone git@github.com:dushoff/$*.git
