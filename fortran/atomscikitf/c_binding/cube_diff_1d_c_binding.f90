module cube_diff_1d_c_binding
    ! Usage:
    use iso_c_binding

    use c_f_string_c_binding

    use askitf_cube_diff_1d, only : cube_diff_1d

    implicit none

    contains 
    
    subroutine c_cube_diff_1d(cube_file_i, cube_file_ii, cube_file_iii) bind(c)


        type(c_ptr), target, intent(in) :: cube_file_i, cube_file_ii, cube_file_iii

        character(len=128), dimension(3) :: cube_file_in_iii
        
        call c_f_string(c_loc(cube_file_i), cube_file_in_iii(1))
        call c_f_string(c_loc(cube_file_ii), cube_file_in_iii(2))
        call c_f_string(c_loc(cube_file_iii), cube_file_in_iii(3))

        call cube_diff_1d(cube_file_in_iii)
    end subroutine c_cube_diff_1d
end module cube_diff_1d_c_binding
