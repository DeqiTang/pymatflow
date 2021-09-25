#
#

CC = gcc
CXX = g++
FC = gfortran #mpif90

CPPFLAGS = 
FPPFLAGS = 
CFLAGS = 
CXXFLAGS = 
FCFLAGS = -g -fopenmp -fPIC

LDFLAGS = 
LDFLAGS += -fopenmp

LIBS +=

modsrc = $(wildcard ./*.f90)
modobj = $(patsubst %.f90, %.o, $(modsrc))

all: $(modobj) lib

lib: $(modobj)
	ar crv libatomsciflowf.a $(modobj) $(LIBS)
	$(FC) -shared -o libatomsciflowf.so $(FCFLAGS) $(FPPFLAGS) $(modobj) $(LIBS) $(LDFLAGS)

%.o: %.f90 $(modobj)
	$(FC) -c $(FPPFLAGS) $(FCFLAGS) $< -o $@

$(modobj): $(modsrc)
	#$(FC) $(FPPFLAGS) $(FCFLAGS) -c $^
	$(shell $(FC) $(FPPFLAGS) $(FCFLAGS) -c *.f90)

install: 
	#
	#
.PHONY: clean

clean:
	rm -rf *.x *.mod *.a *.so *.o $(obj) core 
