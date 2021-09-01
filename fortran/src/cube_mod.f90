module askit_cube_mod
    ! Usage:
    ! askit-cube-1d.x CUBEFILE
    use askit_crystal_mod, only : crystal
    use askit_crystal_mod, only : element_map
    use askit_constant_mod
    
    implicit none

    type :: cube
        type(crystal) :: cube_crystal
        integer :: ngridx, ngridy, ngridz
        real, allocatable :: data(:, :, :)
        contains
        !procedure :: read_cube_file
        procedure, pass :: read_cube_file
    end type
    


    
    contains
    subroutine read_cube_file(this, filename)
        implicit none
        class(cube) :: this
        character(len=*), intent(in) :: filename

        !real, parameter :: bohr_to_angstrom = 0.529177249

        type(element_map) :: ele_map 

        real :: tmp_real

        integer :: i, j, k, tmp_int, last_z_n_data, least_n_six, m
        integer :: ngridx, ngridy, ngridz
        real, allocatable :: data(:)
    
        
        ! read in the input cube file
        open(unit=10, file=filename, status='old', action="read")
    
        do i = 1, 2
            read(10, *)
        end do
    
        read(10, *) this%cube_crystal%natom
        read(10, *) ngridx, this%cube_crystal%cell(1, 1), this%cube_crystal%cell(1, 2), this%cube_crystal%cell(1, 3)
        read(10, *) ngridy, this%cube_crystal%cell(2, 1), this%cube_crystal%cell(2, 2), this%cube_crystal%cell(2, 3)
        read(10, *) ngridz, this%cube_crystal%cell(3, 1), this%cube_crystal%cell(3, 2), this%cube_crystal%cell(3, 3)
    
        this%ngridx = ngridx
        this%ngridy = ngridy
        this%ngridz = ngridz

        do j = 1, 3
            this%cube_crystal%cell(1, j) = this%cube_crystal%cell(1, j) * ngridx * bohr_to_angstrom
        end do
    
        do j = 1, 3
            this%cube_crystal%cell(2, j) = this%cube_crystal%cell(2, j) * ngridy * bohr_to_angstrom
        end do
        
        do j = 1, 3
            this%cube_crystal%cell(3, j) = this%cube_crystal%cell(3, j) * ngridz * bohr_to_angstrom
        end do

        allocate(this%cube_crystal%name(this%cube_crystal%natom))
        allocate(this%cube_crystal%xyz(this%cube_crystal%natom, 3))
        allocate(data(ngridx * ngridy * ngridz))
        allocate(this%data(ngridx, ngridy, ngridz))
    

        call ele_map%initialize()

        do i = 1, this%cube_crystal%natom
            !read(10, *) tmp_int, tmp_real, this%cube_crystal%xyz(i, 1), this%cube_crystal%xyz(i, 2), this%cube_crystal%xyz(i, 3)
            read(10, *) tmp_int, tmp_real, this%cube_crystal%xyz(i, :)
            this%cube_crystal%xyz(i, 1) = this%cube_crystal%xyz(i, 1) * bohr_to_angstrom
            this%cube_crystal%xyz(i, 2) = this%cube_crystal%xyz(i, 2) * bohr_to_angstrom
            this%cube_crystal%xyz(i, 3) = this%cube_crystal%xyz(i, 3) * bohr_to_angstrom
            this%cube_crystal%name(i) = ele_map%get_element_symbol(tmp_int)
            !write(*, *) this%cube_crystal%name(i), this%cube_crystal%xyz(i, :)
        end do
            
        !do i = 1, ngridx * ngridy * ngridz / 6
        !    read(10, *) data((i-1)*6 + 1), data((i-1)*6 + 2), data((i-1)*6 + 3), &
        !        & data((i-1)*6 + 4), data((i-1)*6 + 5), data((i-1)*6 + 6)
        !end do
        ! array in fortran is with column priority, so we cannot reshape 1d data to get 3d data like:
        ! data_3d = reshape(data_1d, (/ngridx, ngridy, ngridz/))
        !do i = 1, ngridx
        !    do j = 1, ngridy
        !        do k = 1, ngridz
        !            this%data(i, j, k) = data((i-1)*ngridy*ngridz + (j-1)*ngridx + k - 1)
        !        end do
        !    end do
        !end do

        least_n_six = ceiling(ngridz / 6.0) - 1
        if (mod(ngridz, 6) == 0) then
            last_z_n_data = 6
        else
            last_z_n_data = mod(ngridz, 6)
        end if

        do i = 1, ngridx
            do j = 1, ngridy
                do k = 1, least_n_six
                    read(10, *) (this%data(i, j, (k-1)*6+m), m=1, 6)
                end do
                read(10, *) (this%data(i, j, least_n_six * 6 + m), m=1, last_z_n_data)
            end do
        end do

        close(unit=10)
    end subroutine

end module askit_cube_mod
