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


all: lib $(modobj) $(obj)

lib: $(obj)
	make -C ../src
	ar crv libbasiccbinding.a $(obj) $(LIBS)

%.o: %.f90
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
