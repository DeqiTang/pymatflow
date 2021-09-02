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

LIBS += ../src/libbasic.a

modsrc = $(wildcard ./*_mod.f90)
modobj = $(patsubst %.f90, %.o, $(modsrc))

src = $(wildcard ./*.f90)
obj = $(patsubst %.f90, %.o, $(src))


all: $(modobj) lib askit-cube-1d.x askit-cube-diff-1d.x askit-cube-to-vtk.x

askit-cube-1d.x: $(obj)
	#ar crv ./src/$(target) $(obj) $(LDFLAGS)
	#$(FC) -o $@	$(FCFLAGS) $(FPPFLAGS) $(obj) $(LIBS) $(LDFLAGS) 
	$(FC) -o $@	$(FCFLAGS) $(FPPFLAGS) cube_1d_program.o $(LIBS) $(LDFLAGS) 

askit-cube-diff-1d.x: $(obj)
	#ar crv ./src/$(target) $(obj) $(LDFLAGS)
	#$(FC) -o $@	$(FCFLAGS) $(FPPFLAGS) $(obj) $(LIBS) $(LDFLAGS) 
	$(FC) -o $@	$(FCFLAGS) $(FPPFLAGS) cube_diff_1d_program.o $(LIBS) $(LDFLAGS) 

askit-cube-to-vtk.x: $(obj)
	$(FC) -o $@ $(FCFLAGS) $(FPPFLAGS) cube_to_vtk_program.o $(LIBS) $(LDFLAGS)

lib: $(obj) src
	make -C ../src
src:
	make -C ../src

%.o: %.f90 src
	$(FC) -c $(FPPFLAGS) $(FCFLAGS) $< -o $@

$(modobj): $(modsrc)
	#$(FC) $(FPPFLAGS) $(FCFLAGS) -c $^
	$(shell $(FC) $(FPPFLAGS) $(FCFLAGS) -c *_mod.f90)

install: 
	cp askit-cube-1d.x ${HOME}/.local/bin
	cp askit-cube-diff-1d.x ${HOME}/.local/bin
	cp askit-cube-to-vtk.x ${HOME}/.local/bin

.PHONY: clean

clean:
	make clean -C ../src
	rm -rf *.x *.mod *.a *.so *.o $(obj) core 
