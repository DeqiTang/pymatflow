module cube_to_vtk_c_binding
    ! Usage:
    ! Reference:
    !   http://fortranwiki.org/fortran/show/c_interface_module
    use iso_c_binding
    use c_f_string_c_binding, only : c_f_string
    use askitf_cube_to_vtk, only : cube_to_vtk

    implicit none

    contains
    subroutine c_cube_to_vtk(cube_file_in_c, vtk_file_out_c) bind(c)
        type(c_ptr), target, intent(in) :: cube_file_in_c, vtk_file_out_c
        !character(c_char), intent(in) :: cube_file_in_c, vtk_file_out_c
        character(len=128) :: cube_file_in, vtk_file_out

        call c_f_string(c_loc(cube_file_in_c), cube_file_in)
        call c_f_string(c_loc(vtk_file_out_c), vtk_file_out)

        call cube_to_vtk(cube_file_in, vtk_file_out)
    end subroutine c_cube_to_vtk
end module cube_to_vtk_c_binding
