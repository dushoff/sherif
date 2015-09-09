### Hooks for the editor to set the default target
current: target

target pngtarget pdftarget vtarget acrtarget: sim.exe 

##################################################################

# make files

Sources = Makefile .gitignore

##################################################################

# Compiler
# CC = g++
# CC_FLAGS = -Wall -O3 -std=c++11 -Wno-sign-compare

CC := clang++
CC_FLAGS := -Wall -O3 -Wno-predefined-identifier-outside-function

code = mc.cpp simulator.cpp individual.cpp dcTools.cpp RV.cpp globalVar.cpp dcMatrix.cpp EventNumberContainer.cpp 
# Sources += $(code)

obj = $(code:%.cpp=%.o)
header = $(code:%.cpp=%.h)

$(obj): %.o: %.cpp %.h
	$(CC) $(CC_FLAGS) -c $*.cpp

# This is a tangle, and I am importing all of the .h files for now.
# I think it's better style to _not_ include .h files in .h files (each .cpp should include what it needs directly).
main.o: individual.h mc.h
mc.o: globalVar.h simulator.h
RV.o: dcTools.h globalVar.h
individual.o: simulator.h EventNumberContainer.h dcTools.h
simulator.o: individual.h

sim.exe: sim.cpp
%.exe: %.o $(obj)
	$(CC) $(CC_FLAGS) $(MAIN) $(CPP) -o $@

### Temp ###
$(code):
	/bin/cp ~/Dropbox/SHERIF/$@ .

$(header):
	/bin/cp ~/Dropbox/SHERIF/$@ .

sim.cpp:
	/bin/cp ~/Dropbox/SHERIF/main.cpp .

Makefile: $(header)

######################################################################

# JD git rules

ms = ../makestuff
include $(ms)/git.mk

