#
#

CC = gcc
CXX = g++
FC = gfortran #mpif90

CPPFLAGS = 
FPPFLAGS = -I ../src
CFLAGS = 
CXXFLAGS = 
FCFLAGS = -g -fopenmp -fPIC

LDFLAGS = 
LDFLAGS += -fopenmp

LIBS += ../src/libatomsciflowf.a

src = $(wildcard ./*.f90)
obj = $(patsubst %.f90, %.o, $(src))


all: lib asflowf-cube-1d.x asflowf-cube-diff-1d.x asflowf-cube-to-vtk.x

asflowf-cube-1d.x: $(obj)
	#ar crv ./src/$(target) $(obj) $(LDFLAGS)
	#$(FC) -o $@	$(FCFLAGS) $(FPPFLAGS) $(obj) $(LIBS) $(LDFLAGS) 
	$(FC) -o $@	$(FCFLAGS) $(FPPFLAGS) cube_1d_program.o $(LIBS) $(LDFLAGS) 

asflowf-cube-diff-1d.x: $(obj)
	#ar crv ./src/$(target) $(obj) $(LDFLAGS)
	#$(FC) -o $@	$(FCFLAGS) $(FPPFLAGS) $(obj) $(LIBS) $(LDFLAGS) 
	$(FC) -o $@	$(FCFLAGS) $(FPPFLAGS) cube_diff_1d_program.o $(LIBS) $(LDFLAGS) 

asflowf-cube-to-vtk.x: $(obj)
	$(FC) -o $@ $(FCFLAGS) $(FPPFLAGS) cube_to_vtk_program.o $(LIBS) $(LDFLAGS)

lib: $(obj) src
	make -C ../src
src:
	make -C ../src

%.o: %.f90 src
	$(FC) -c $(FPPFLAGS) $(FCFLAGS) $< -o $@


install: 
	cp asflowf-cube-1d.x ${HOME}/.local/bin
	cp asflowf-cube-diff-1d.x ${HOME}/.local/bin
	cp asflowf-cube-to-vtk.x ${HOME}/.local/bin

.PHONY: clean

clean:
	make clean -C ../src
	rm -rf *.x *.mod *.a *.so *.o $(obj) core 
