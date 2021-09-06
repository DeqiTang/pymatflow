module cube_1d_c_binding
    ! Usage:
    use iso_c_binding

    use c_f_string_c_binding, only : c_f_string

    use askitf_cube_1d, only : cube_1d

    implicit none
    
    !real, parameter :: bohr_to_angstrom = 0.529177249

    contains 
    
    subroutine c_cube_1d(cube_file_in) bind(c)
        type(c_ptr), target, intent(in) :: cube_file_in
        character(len=128) :: cube_file_in_f
        
        call c_f_string(c_loc(cube_file_in), cube_file_in_f)

        call cube_1d(cube_file_in_f)
    end subroutine c_cube_1d
end module cube_1d_c_binding
