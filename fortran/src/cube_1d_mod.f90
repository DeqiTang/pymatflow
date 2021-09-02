module cube_1d_mod

    use askit_crystal_mod
    use askit_cube_mod
    use askit_constant_mod
    
    implicit none

    contains

    subroutine cube_1d(cube_file_in)
        integer :: i, j, k
        integer :: ngridx, ngridy, ngridz

        real, allocatable :: data_red_a(:), data_red_b(:), data_red_c(:)
        type(cube) :: cube_i

        real :: cell_volume, cell_volume_per_unit, tmp, tmp_vec(3)
        real :: a, b, c
        real :: total_electron

        ! character, allocatable :: cube_file_in
        character(len=128), intent(in) :: cube_file_in

        ! command line output 
        write(*, *) "*******************************************************************************"
        write(*, *) "***                 CUBE FILE PROCESSOR FROM ATOMSCIKIT                     ***"
        write(*, *) "*******************************************************************************"


        ! read cube file
        write(*, *) "On getting the command line argument:"
        if ( cube_file_in == "") then
            write(*, *) "You should provide the name for the input cube file!"
            stop
        else
            write(*, *) "The input cube file name is: ", cube_file_in
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

        allocate(data_red_a(ngridx))
        allocate(data_red_b(ngridy))
        allocate(data_red_c(ngridz))

        !$OMP PARALLEL
        !$OMP DO
        do i = 1, ngridx
            data_red_a(i) = 0
            do j = 1, ngridy
                do k = 1, ngridz
                    data_red_a(i) =  data_red_a(i) + cube_i%data(i, j, k)
                end do
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        data_red_a = data_red_a * cell_volume_per_unit / bohr_to_angstrom**3

        !$OMP PARALLEL
        !$OMP DO
        do j = 1, ngridy
            data_red_b(j) = 0
            do i = 1, ngridx
                do k = 1, ngridz
                    data_red_b(j) = data_red_b(j) + cube_i%data(i, j, k)
                end do
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        data_red_b = data_red_b * cell_volume_per_unit / bohr_to_angstrom**3

        !$OMP PARALLEL
        !$OMP DO
        do k = 1, ngridz
            data_red_c(k) = 0
            do i = 1, ngridx
                do j = 1, ngridy
                    data_red_c(k) = data_red_c(k) + cube_i%data(i, j, k)
                end do
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        data_red_c = data_red_c * cell_volume_per_unit / bohr_to_angstrom**3

        call dot_3(cube_i%cube_crystal%cell(1, :), cube_i%cube_crystal%cell(1, :), a)
        a = sqrt(a)
        call dot_3(cube_i%cube_crystal%cell(2, :), cube_i%cube_crystal%cell(2, :), b)
        b = sqrt(b)
        call dot_3(cube_i%cube_crystal%cell(3, :), cube_i%cube_crystal%cell(3, :), c)
        c = sqrt(c)


        ! output dimension reduced data
        open(unit=11, file="charge.1d.c.data", action="write")
        write(11, *) "#c(angstrom) rho(e) (number of electron per Angstrom)"
        do i = 1, ngridz
            write(11, *) (0 + c / (ngridz - 1) * (i - 1)), data_red_c(i)
        end do
        close(unit=11)
        
        open(unit=11, file="charge.1d.b.data", action="write")
        write(11, *) "#b(angstrom) rho(e) (number of electron per Angstrom)"
        do i = 1, ngridy
            write(11, *) (0 + b / (ngridy - 1) * (i - 1)), data_red_b(i)
        end do
        close(unit=11)

        open(unit=11, file="charge.1d.a.data", action="write")
        write(11, *) "#a(angstrom) rho(e) (number of electron per Angstrom)"
        do i = 1, ngridx
            write(11, *) (0 + a / (ngridx - 1) * (i - 1)), data_red_a(i)
        end do
        close(unit=11)

        ! output the total structure
        call write_xyz(cube_i%cube_crystal, "charge-1d-structure.xyz")


        stop
    end subroutine cube_1d

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
end module cube_1d_mod
