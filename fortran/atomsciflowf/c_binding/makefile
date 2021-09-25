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

modsrc = $(wildcard ./*.f90)
modobj = $(patsubst %.f90, %.o, $(modsrc))

all: lib $(modobj)

lib:
	make -C ../src
	ar crv libatomsciflowf-c-binding.a $(obj) $(LIBS)

%.o: %.f90
	$(FC) -c $(FPPFLAGS) $(FCFLAGS) $< -o $@

$(modobj): $(modsrc)
	#$(FC) $(FPPFLAGS) $(FCFLAGS) -c $^
	$(shell $(FC) $(FPPFLAGS) $(FCFLAGS) -c *.f90)

install: 
	cp asflow-cube-1d.x ${HOME}/.local/bin
	cp asflow-cube-diff-1d.x ${HOME}/.local/bin
	cp asflow-cube-to-vtk.x ${HOME}/.local/bin

.PHONY: clean

clean:
	make clean -C ../src
	rm -rf *.x *.mod *.a *.so *.o $(obj) core 
