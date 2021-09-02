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

modsrc = $(wildcard ./*_mod.f90)
modobj = $(patsubst %.f90, %.o, $(modsrc))

src = $(wildcard ./*.f90)
obj = $(patsubst %.f90, %.o, $(src))


all: $(modobj) lib

lib: $(obj)
	ar crv libaskitf.a $(obj) $(LIBS)
	$(FC) -shared -o libaskitf.so $(FCFLAGS) $(FPPFLAGS) $(obj) $(LIBS) $(LDFLAGS)

%.o: %.f90 $(modobj)
	$(FC) -c $(FPPFLAGS) $(FCFLAGS) $< -o $@

$(modobj): $(modsrc)
	#$(FC) $(FPPFLAGS) $(FCFLAGS) -c $^
	$(shell $(FC) $(FPPFLAGS) $(FCFLAGS) -c *_mod.f90)

install: 
	#
	#
.PHONY: clean

clean:
	rm -rf *.x *.mod *.a *.so *.o $(obj) core 
