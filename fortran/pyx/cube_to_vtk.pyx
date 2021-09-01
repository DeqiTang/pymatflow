#cdef extern from "pygfunc.h":
#    void c_gfunc(double* a, int* n, int* m, double* a, double* b, double* c)

cdef extern:
    void c_cube_to_vtk(char cube_file_in, char cube_file_out)


def cube_to_vtk(char cube_file_in, char vtk_file_out):
    cdef:
        ar[double] ax = linspace(a, b, n)
        ar[double,ndim=2] c = empty((n, n), order='F')
    c_cube_to_vtk(&cube_file_in, &vtk_file_out)