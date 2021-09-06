module askitf_cube_diff_1d

    use askitf_crystal, only : write_xyz
    use askitf_cube, only : cube
    use askitf_constants

    implicit none
    

    contains 
    
    subroutine cube_diff_1d(cube_file_in_iii)

        type(cube), dimension(3) :: cube_iii

        integer :: i, j, k
        integer :: natom, ngridx, ngridy, ngridz
        real(kind=dp), allocatable :: data_diff(:, :, :), data_red_a(:), data_red_b(:), data_red_c(:)

        real(kind=dp) :: cell_volume, cell_volume_per_unit, tmp, tmp_vec(3)
        real(kind=dp) :: a, b, c
        real(kind=dp) :: total_electron

        character(len=128), dimension(3), intent(in) :: cube_file_in_iii

        ! command line output 
        write(*, *) "*******************************************************************************"
        write(*, *) "***                 CUBE FILE PROCESSOR FROM ATOMSCIKIT                     ***"
        write(*, *) "*******************************************************************************"


        ! read cube file 

        if ( cube_file_in_iii(1) == "") then
            write(*, *) "You should provide the name for the input cube file!"
            stop
        else
            write(*, *) "The input cube file name is: "
            write(*, *) cube_file_in_iii(1)
            write(*, *) cube_file_in_iii(2)
            write(*, *) cube_file_in_iii(3)
        end if

        ! read in the input cube file
        call cube_iii(1)%read_cube_file(cube_file_in_iii(1))
        call cube_iii(2)%read_cube_file(cube_file_in_iii(2))
        call cube_iii(3)%read_cube_file(cube_file_in_iii(3))

        ngridx = cube_iii(1)%ngridx
        ngridy = cube_iii(1)%ngridy
        ngridz = cube_iii(1)%ngridz

        ! value in cube file are \rho(r)_of_electrons in unit of e/Bohr^3
        ! namely number of electrons each Borh^3
        ! so we have to convert it to e/Angstrom^3, through divide it by bohr_to_angstrom**3

        call cross_3(cube_iii(1)%cube_crystal%cell(1, :), cube_iii(1)%cube_crystal%cell(2, :), tmp_vec)

        call dot_3(tmp_vec, cube_iii(1)%cube_crystal%cell(3, :), cell_volume)

        cell_volume_per_unit = cell_volume / ngridx / ngridy / ngridz

        data_diff = cube_iii(1)%data - cube_iii(2)%data - cube_iii(3)%data

        total_electron = sum(cube_iii(1)%data) * cell_volume_per_unit / bohr_to_angstrom**3

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

        !$OMP PARALLEL DO
        do i = 1, ngridx
            data_red_a(i) = 0
            do j = 1, ngridy
                do k = 1, ngridz
                    data_red_a(i) =  data_red_a(i) + data_diff(i, j, k)
                end do
            end do
        end do
        !$OMP END PARALLEL DO
        data_red_a = data_red_a * cell_volume_per_unit / bohr_to_angstrom**3

        !$OMP PARALLEL
        !$OMP DO
        do j = 1, ngridy
            data_red_b(j) = 0
            do i = 1, ngridx
                do k = 1, ngridz
                    data_red_b(j) = data_red_b(j) + data_diff(i, j, k)
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
                    data_red_c(k) = data_red_c(k) + data_diff(i, j, k)
                end do
            end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        data_red_c = data_red_c * cell_volume_per_unit / bohr_to_angstrom**3

        call dot_3(cube_iii(1)%cube_crystal%cell(1, :), cube_iii(1)%cube_crystal%cell(1, :), a)
        a = sqrt(a)
        call dot_3(cube_iii(1)%cube_crystal%cell(2, :), cube_iii(1)%cube_crystal%cell(2, :), b)
        b = sqrt(b)
        call dot_3(cube_iii(1)%cube_crystal%cell(3, :), cube_iii(1)%cube_crystal%cell(3, :), c)
        c = sqrt(c)


        ! output dimension reduced data
        open(unit=11, file="diff.charge.1d.c.data", action="write")
        write(11, *) "#c(angstrom) \delta(rho(e)) (number of electron per Angstrom)"
        do i = 1, ngridz
            write(11, *) (0 + c / (ngridz - 1) * (i - 1)), data_red_c(i)
        end do
        close(unit=11)
        
        open(unit=11, file="diff.charge.1d.b.data", action="write")
        write(11, *) "#b(angstrom) \delta(rho(e)) (number of electron per Angstrom)"
        do i = 1, ngridy
            write(11, *) (0 + b / (ngridy - 1) * (i - 1)), data_red_b(i)
        end do
        close(unit=11)

        open(unit=11, file="diff.charge.1d.a.data", action="write")
        write(11, *) "#a(angstrom) \delta(rho(e)) (number of electron per Angstrom)"
        do i = 1, ngridx
            write(11, *) (0 + a / (ngridx - 1) * (i - 1)), data_red_a(i)
        end do
        close(unit=11)

        ! output the total structure
        call write_xyz(cube_iii(1)%cube_crystal, "diff-charge-total.xyz")
        
        stop
    end subroutine cube_diff_1d

    subroutine cross_3(x, y, z)
        implicit none
        real(kind=dp), dimension(3), intent(in) :: x, y
        real(kind=dp), dimension(3), intent(out) :: z

        z(1) = x(2) * y(3) - x(3) * y(2)
        z(2) = x(3) * y(1) - x(1) * y(3)
        z(3) = x(1) * y(2) - x(2) * y(1)
    end subroutine

    subroutine dot_3(x, y, z)
        implicit none
        real(kind=dp), dimension(3), intent(in) :: x, y 
        real(kind=dp), intent(out) :: z 
        z = x(1) * y(1) + x(2) * y(2) + x(3) * y(3)
    end subroutine dot_3
end module askitf_cube_diff_1d
