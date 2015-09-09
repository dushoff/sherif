### Hooks for the editor to set the default target
current: target

target pngtarget pdftarget vtarget acrtarget: sim.exe 

##################################################################

# make files

Sources = Makefile .gitignore 

##################################################################

# compiler.mk is intended to be a local file; it uses clang by default, but you can just edit and save it on any particular machine. Don't add it to sources

compiler.mk:
	/bin/cp clang.mk $@
Sources += clang.mk gcc.mk

# Old C stuff that we may not need
Sources += exe.mk

### Temp ###
$(code):
	/bin/cp ~/Dropbox/SHERIF/$@ .

$(header):
	/bin/cp ~/Dropbox/SHERIF/$@ .

######################################################################

# JD git rules

ms = ../makestuff
-include $(ms)/git.mk
