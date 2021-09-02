module cube_to_vtk_mod

    use askit_crystal_mod
    use askit_cube_mod
    use askit_constant_mod
    
    implicit none

    contains 

    subroutine cube_to_vtk(cube_file_in, vtk_file_out)
            
        integer :: i, j, k
        integer :: ngridx, ngridy, ngridz

        type(cube) :: cube_i

        real :: cell_volume, cell_volume_per_unit, tmp, tmp_vec(3)
        real :: a, b, c, x, y, z, total_electron

        ! character, allocatable :: cube_file_in
        character(len=128), intent(in) :: cube_file_in, vtk_file_out

        ! command line output 
        write(*, *) "*******************************************************************************"
        write(*, *) "***                 CUBE FILE PROCESSOR FROM ATOMSCIKIT                     ***"
        write(*, *) "*******************************************************************************"


        ! read cube file
        
        write(*, *) "On getting the command line argument:"
        if ( cube_file_in == "" .or. vtk_file_out == "") then
            write(*, *) "You should provide the name for the input cube file and output vtk file!"
            stop
        else
            write(*, *) "The input cube file name is: ", cube_file_in
            write(*, *) "The output vtk file name is: ", vtk_file_out
        end if

        call cube_i%read_cube_file(cube_file_in)

        write(*, *) "Successfully read the cube file!"

        ! value in cube file are \rho(r)_of_electrons in unit of e/Bohr^3
        ! namely number of electrons each Borh^3
        ! so we have to convert it to e/Angstrom^3, through divide it by bohr_to_angstrom**3

        call cross_3(cube_i%cube_crystal%cell(1, :), cube_i%cube_crystal%cell(2, :), tmp_vec)

        call dot_3(tmp_vec, cube_i%cube_crystal%cell(3, :), cell_volume)

        ngridx = cube_i%ngridx
        ngridy = cube_i%ngridy
        ngridz = cube_i%ngridz

        cell_volume_per_unit = cell_volume / ngridx / ngridy / ngridz

        total_electron = sum(cube_i%data) * cell_volume_per_unit / bohr_to_angstrom**3
        
        write(*, *) "-----------------------------------------------------------------"
        write(*, *) "                 Out put collected information                   "
        write(*, *) "-----------------------------------------------------------------"
        write(*, *) "ngridx: ", ngridx
        write(*, *) "ngridy: ", ngridy
        write(*, *) "ngridz: ", ngridz
        write(*, *) "cell volume: ", cell_volume
        write(*, *) "total number of electrons: ", total_electron
        write(*, *) "-----------------------------------------------------------------"

        call dot_3(cube_i%cube_crystal%cell(1, :), cube_i%cube_crystal%cell(1, :), a)
        a = sqrt(a)
        call dot_3(cube_i%cube_crystal%cell(2, :), cube_i%cube_crystal%cell(2, :), b)
        b = sqrt(b)
        call dot_3(cube_i%cube_crystal%cell(3, :), cube_i%cube_crystal%cell(3, :), c)
        c = sqrt(c)

        !write(*, *) a, b ,c
        
        ! output data to vtk file
        open(11, file=vtk_file_out, status="replace", action="write")
        write(11, "(A)") "# vtk DataFile Version 5.1"
        write(11, "(A)") cube_file_in
        write(11, "(A)") "ASCII"
        write(11, "(A)") "DATASET STRUCTURED_POINTS"
        write(11, "(A, 3I10)") "DIMENSIONS ", ngridx, ngridy, ngridz
        write(11, "(A, 3F15.6)") "SPACING", a/ngridx, b/ngridy, c/ngridz


        !write(11, "(A, I10, A)") "POINTS ", ngridx * ngridy * ngridz, ' float'
        !do i = 1, ngridx
        !    do j = 1, ngridy
        !        do k = 1, ngridz
        !            !x = real(i-1) / real(ngridx) * a
        !            !y = real(j-1) / real(ngridy) * b
        !            !z = real(k-1) / real(ngridz) * c
        !            x = real(i) / real(ngridx) * a
        !            y = real(j) / real(ngridy) * b
        !            z = real(k) / real(ngridz) * c
        !            write(11, *) x, y, z
        !        end do
        !    end do
        !end do
        
        ! write grid point values:
        write(11, "(A, I10)") 'POINT_DATA ', ngridx * ngridy * ngridz
        write(11, "(A)") 'SCALARS CON float 1'
        write(11, "(A)") 'LOOKUP_TABLE default'
        do k = 1, ngridz
            do j = 1, ngridy
                do i = 1, ngridx
                    write(11, *) cube_i%data(i, j, k)
                end do
            end do
        end do
        close(11)
        ! output the total structure
        call write_xyz(cube_i%cube_crystal, "cube-structure.xyz")
    end subroutine cube_to_vtk


    subroutine cross_3(x, y, z)
        implicit none
        real, dimension(3), intent(in) :: x, y
        real, dimension(3), intent(out) :: z

        z(1) = x(2) * y(3) - x(3) * y(2)
        z(2) = x(3) * y(1) - x(1) * y(3)
        z(3) = x(1) * y(2) - x(2) * y(1)
    end subroutine

    subroutine dot_3(x, y, z)
        implicit none
        real, dimension(3), intent(in) :: x, y 
        real, intent(out) :: z 
        z = x(1) * y(1) + x(2) * y(2) + x(3) * y(3)
    end subroutine dot_3
end module cube_to_vtk_mod
