program cube_1d_program
    ! Usage:
    ! askit-cube-1d.x CUBEFILE

    use asflowf_cube_1d, only : cube_1d

    implicit none

    ! character, allocatable :: cube_file_in
    character(len=128) :: cube_file_in
    
    ! read cube file specified by the first command argument
    call get_command_argument(1, cube_file_in)

    call cube_1d(cube_file_in)
    
    stop
end program cube_1d_program
