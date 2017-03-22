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
	./gridgen 1
	$(root)/bin/solver
	mv mpnt.dat mpnt_10.dat
	mv data_out.cgns data_out_10.cgns
	./gridgen 2
	$(root)/bin/solver
	mv mpnt.dat mpnt_20.dat
	mv data_out.cgns data_out_20.cgns
	./gridgen 4
	$(root)/bin/solver
	mv mpnt.dat mpnt_40.dat
	mv data_out.cgns data_out_40.cgns
	