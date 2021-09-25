program cube_diff_1d_program
    ! Usage:
    ! askit-cube-diff-1d.x CUBEFILE1 CUBEFILE2 CUBEFILE3
    use asflowf_cube_diff_1d, only : cube_diff_1d

    implicit none
    
    character(len=128), dimension(3) :: cube_file_in_iii

    ! command line output 
    write(*, *) "*******************************************************************************"
    write(*, *) "***                 CUBE FILE PROCESSOR FROM ATOMSCIKIT                     ***"
    write(*, *) "*******************************************************************************"


    ! read cube file specified by the first command argument
    call get_command_argument(1, cube_file_in_iii(1))
    call get_command_argument(2, cube_file_in_iii(2))
    call get_command_argument(3, cube_file_in_iii(3))

    !
    call cube_diff_1d(cube_file_in_iii)

    stop
end program cube_diff_1d_program