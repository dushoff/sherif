### Hooks for the editor to set the default target
current: target

target pngtarget pdftarget vtarget acrtarget: sherif 

##################################################################

# make files

Sources = Makefile 

##################################################################

# C++ stuff

CC = g++
CC_FLAGS = -Wall -O3 -std=c++11 -Wno-sign-compare

CPP = mc.cpp simulator.cpp individual.cpp dcTools.cpp RV.cpp globalVar.cpp dcMatrix.cpp EventNumberContainer.cpp 

Sources += $(CPP) $(MAIN)

sherif: $(MAIN) $(CPP)
	$(CC) $(CC_FLAGS) $(MAIN) $(CPP) -o $@

### Temp ###
$(CPP) $(MAIN):
	/bin/cp ~/Dropbox/SHERIF/$@ .

######################################################################

# JD weird generic stuff

md = ../make/
rrd = ../RR/

local = default
-include $(md)/local.mk
-include $(md)/$(local).mk
-include $(rrd)/inc.mk

