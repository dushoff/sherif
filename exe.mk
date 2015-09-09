
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
%.exe: %.cpp $(code)
	$(CC) $(CC_FLAGS) $^ -o $@


sim.cpp:
	/bin/cp ~/Dropbox/SHERIF/main.cpp $@

Makefile: $(header)
