root=../../../
local=`pwd`
CGNSlibs = /usr/lib/libcgns.so -lhdf5
LIBS = $(CGNSlibs) -lm
all: compile_solver gridgen data_in.cgns
	$(root)/bin/solver
data_in.cgns: gridgen
	./$<


gridgen: gridgen.F90
	gfortran -J$(root)obj $< -o $@ $(LIBS)
#	@cd $(root)/obj
#	gfortran -o $(local)/$@ $(local)$< $(LIBS)

compile_solver:
	$(MAKE) -C $(root) all
	
.PHONY: clean

clean:
	rm -vf *.cgns *.o *.mod *.bin gridgen *.dat
	
	
test: compile_solver gridgen
	make -C $(root)/tools
	./gridgen 100
	$(root)/bin/solver
	$(root)/bin/cgns2line data_out.cgns line_100.dat
	mv data_out.cgns data_out_100.cgns
	./gridgen 50
	$(root)/bin/solver
	$(root)/bin/cgns2line data_out.cgns line_50.dat
	mv data_out.cgns data_out_50.cgns
	./gridgen 10
	$(root)/bin/solver
	$(root)/bin/cgns2line data_out.cgns line_10.dat
	mv data_out.cgns data_out_10.cgns
	