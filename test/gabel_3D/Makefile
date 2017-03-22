root=../../
local=`pwd`
CGNSlibs = /usr/lib/libcgns.so -lhdf5
LIBS = $(CGNSlibs) -lm
all: compile_solver gridgen git.cgns
	$(root)/bin/solver
git.cgns: gridgen
	./$<


gridgen: gridgen.F90
	gfortran -J$(root)obj $< -o $@ $(LIBS)
#	@cd $(root)/obj
#	gfortran -o $(local)/$@ $(local)$< $(LIBS)

compile_solver:
	$(MAKE) -C $(root) all
	
.PHONY: clean

clean:
	rm -vf *.cgns *.o *.mod *.bin gridgen