CC = gcc
CXX = g++-10 #mpicxx
#NVCC = nvcc

#CPPFLAGS = -g 
#CPPFLAGS += -I../include

#CFLAGS = -g
#CXXFLAGS = -g 
#NVCCFLAGS = -g -arch sm_50 ## here use "nvcc -h" to find an allowed or supported value in --gpu-architecture.

#LDFLAGS = -lcppunit -ldl -lcudart -lcufft -lmpi

#src = $(wildcard ./*.cpp ../src/*.cpp)
#obj = $(patsubst %.cpp, %.o, $(src))

#src_cu = $(wildcard ./*.cu ../src/*.cu)
#obj_cu = $(patsubst %.cu, %.cu.o, $(src_cu))

#target = main.exe

#$(target): $(obj) $(obj_cu)
#	$(CXX) -o $(target) $(CPPFLAGS) $(obj) $(obj_cu) $(LDFLAGS)

#%.o: %.cpp
#	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@


#%.cu.o: %.cu
#	$(NVCC) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@


#.PHONY: clean

#clean:
#	rm -rf $(obj) $(obj_cu) $(target) core

get_pes:
	$(CXX) -o get_pes.x get_pes.cpp -std=c++17

install: get_pes
	cp get_pes.x ${HOME}/.local/bin