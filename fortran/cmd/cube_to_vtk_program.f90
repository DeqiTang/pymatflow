program cube_to_vtk_program
    ! Usage:
    ! askit-cube-to-vtk.x CUBEFILE VTKFILE
    use cube_to_vtk_mod, only : cube_to_vtk

    implicit none

    character(len=128) :: cube_file_in, vtk_file_out

    ! read cube file specified by the first command argument
    call get_command_argument(1, cube_file_in)
    call get_command_argument(2, vtk_file_out)

    
    ! 
    call cube_to_vtk(cube_file_in, vtk_file_out)

    stop
end program cube_to_vtk_program
