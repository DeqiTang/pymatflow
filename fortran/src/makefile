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


all: $(modobj) askit-cube-1d.x askit-cube-diff-1d.x askit-cube-to-vtk.x lib

askit-cube-1d.x: $(obj)
	#ar crv ./src/$(target) $(obj) $(LDFLAGS)
	#$(FC) -o $@	$(FCFLAGS) $(FPPFLAGS) $(obj) $(LIBS) $(LDFLAGS) 
	$(FC) -o $@	$(FCFLAGS) $(FPPFLAGS) crystal_mod.o cube_mod.o cube_1d.o $(LIBS) $(LDFLAGS) 

askit-cube-diff-1d.x: $(obj)
	#ar crv ./src/$(target) $(obj) $(LDFLAGS)
	#$(FC) -o $@	$(FCFLAGS) $(FPPFLAGS) $(obj) $(LIBS) $(LDFLAGS) 
	$(FC) -o $@	$(FCFLAGS) $(FPPFLAGS) crystal_mod.o cube_mod.o cube_diff_1d.o $(LIBS) $(LDFLAGS) 

askit-cube-to-vtk.x: $(obj)
	$(FC) -o $@ $(FCFLAGS) $(FPPFLAGS) crystal_mod.o cube_mod.o cube_to_vtk.o $(LIBS) $(LDFLAGS)

lib: $(obj)
	ar crv libbasic.a cube_to_vtk_c_binding_mod.o crystal_mod.o cube_mod.o $(LIBS)
	$(FC) -shared -o libbasic.so $(FCFLAGS) $(FPPFLAGS) cube_to_vtk_c_binding_mod.o crystal_mod.o cube_mod.o $(LIBS) $(LDFLAGS)

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
	rm -rf *.x *.mod *.a *.so *.o $(obj) core 
