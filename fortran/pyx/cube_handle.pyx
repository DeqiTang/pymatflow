"""
Reference:
    https://cython-docs2.readthedocs.io/en/latest/src/tutorial/strings.html
"""
#cimport cython


cdef extern:
    void c_cube_to_vtk(char* cube_file_in, char* cube_file_out)

#def cube_to_vtk(char cube_file_in, char vtk_file_out):
def cube_to_vtk(str cube_file_in, str vtk_file_out):
    # this will not compile !
    # cdef char* cube_in = cube_file_in.encode("UTF-8")
    # becuase this takes the pointer to the byte buffer of the Python byte string. 
    # Trying to do the same without keeping a reference to the Python byte string 
    # will fail with a compile error
    py_byte_string = cube_file_in.encode("UTF-8")
    cdef char* cube_in = py_byte_string
    py_byte_string = vtk_file_out.encode("UTF-8")
    cdef char* vtk_out = py_byte_string

    c_cube_to_vtk(cube_in, vtk_out)


cdef extern:
    void c_cube_diff_1d(char* cube_file_i, char* cube_file_ii, char* cube_file_iii)

def cube_diff_1d(str cube_file_i, str cube_file_ii, str cube_file_iii):
    py_byte_string = cube_file_i.encode("UTF-8")
    cdef char* cube_in_i = py_byte_string
    py_byte_string = cube_file_ii.encode("UTF-8")
    cdef char* cube_in_ii = py_byte_string
    py_byte_string = cube_file_iii.encode("UTF-8")
    cdef char* cube_in_iii = py_byte_string

    c_cube_diff_1d(cube_in_i, cube_in_ii, cube_in_iii)

cdef extern:
    void c_cube_1d(char* cube_file_in)

def cube_1d(str cube_file_in):
    py_byte_string = cube_file_in.encode("UTF-8")
    cdef char* cube_in = py_byte_string

    c_cube_1d(cube_in)