module cube_to_vtk_c_binding_mod
    ! Usage:
    ! Reference:
    !   http://fortranwiki.org/fortran/show/c_interface_module
    use iso_c_binding
    use askit_crystal_mod
    use askit_cube_mod
    use askit_constant_mod

    implicit none

    character(len=1,kind=C_char), parameter :: NUL = C_NULL_char

    interface c_f_string
        module procedure c_f_string_chars
        module procedure c_f_string_ptr
    end interface

    contains
    subroutine c_cube_to_vtk(cube_file_in_c, vtk_file_out_c) bind(c)
        integer :: i, j, k
        integer :: ngridx, ngridy, ngridz

        type(cube) :: cube_i

        real :: cell_volume, cell_volume_per_unit, tmp, tmp_vec(3)
        real :: a, b, c, x, y, z, total_electron

        type(c_ptr), target, intent(in) :: cube_file_in_c, vtk_file_out_c
        !character(c_char), intent(in) :: cube_file_in_c, vtk_file_out_c
        character(len=128) :: cube_file_in, vtk_file_out

        call c_f_string(c_loc(cube_file_in_c), cube_file_in)
        call c_f_string(c_loc(vtk_file_out_c), vtk_file_out)

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
    end subroutine c_cube_to_vtk

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
    !
    subroutine c_f_string_chars(c_string, f_string)
        character(len=1, kind=c_char), intent(in) :: c_string(*)
        character(len=*), intent(out) :: f_string
        integer :: i
        i = 1
        do while (c_string(i) /= nul .and. i <= len(f_string))
            f_string(i:i) = c_string(i)
            i = i + 1
        end do
        if (i < len(f_string)) f_string(i:) = ' '
    end subroutine c_f_string_chars

    subroutine c_f_string_ptr(c_string, f_string)
        type(C_ptr), intent(in) :: C_string
        character(len=*), intent(out) :: F_string
        character(len=1,kind=C_char), dimension(:), pointer :: p_chars
        integer :: i
        if (.not. C_associated(C_string)) then
            F_string = ' '

        else
            call C_F_pointer(C_string,p_chars,[huge(0)])
            i=1
            do while(p_chars(i)/=NUL .and. i<=len(F_string))
                F_string(i:i) = p_chars(i)
                i=i+1
            end do
            if (i<len(F_string)) F_string(i:) = ' '
        end if
    end subroutine C_F_string_ptr
end module cube_to_vtk_c_binding_mod
