module askitf_crystal

    use askitf_constants
    use askitf_kinds, only : dp

    implicit none
    
    !type :: atom
    !    real, dimension(3) :: xyz
    !    character(len=128) :: name
    !end type atom
    
    !real, parameter :: bohr_to_angstrom = 0.529177249

    type :: crystal
        real(kind=dp), dimension(3, 3) :: cell
        integer :: natom
        character(len=2), dimension(:), allocatable :: name
        real(kind=dp), dimension(:, :), allocatable :: xyz
    end type crystal

    type :: element
        integer :: number
        real(kind=dp) :: mass
        character(len=2) :: symbol
    end type element

    type :: element_map
        type(element), dimension(118) :: elements
        contains
        procedure :: initialize => initialize_element_map
        procedure :: get_element_symbol
        procedure :: get_element_number
    end type

    contains
    subroutine read_xyz(cryst, filename)
        implicit none

        type(crystal), intent(inout) :: cryst
        character(len=*) :: filename

        type(element_map) :: ele_map

        integer :: i, j, tmp_int
        character(len=128) :: tmp_str_vec(4)

        open(unit=101, file=filename, status="OLD", action="READ")

        read(101, *) cryst%natom
        read(101, *) tmp_str_vec(1), &
            & cryst%cell(1, 1), cryst%cell(1, 2), cryst%cell(1, 3), &
            & tmp_str_vec(2), &
            & cryst%cell(2, 1), cryst%cell(2, 2), cryst%cell(2, 3),&
            & tmp_str_vec(3), &
            &cryst%cell(3, 1), cryst%cell(3, 2), cryst%cell(3, 3)
        
        allocate(cryst%xyz(cryst%natom, 3))
        
        call ele_map%initialize()

        do i = 1, cryst%natom
            read(101, *) tmp_int, cryst%xyz(i, :)
            cryst%name(i) = ele_map%get_element_symbol(tmp_int)
        end do
        close(unit=101)
    end subroutine

    subroutine write_xyz(cryst, filename)
        implicit none

        type(crystal), intent(in) :: cryst
        character(len=*), intent(in) :: filename

        integer :: i, j, k

        open(unit=102, file=filename, status="UNKNOWN", action="WRITE")
        write(102, *) cryst%natom
        write(102, *) "cell: ", cryst%cell(1, 1), " ", cryst%cell(1, 2), " ", cryst%cell(1, 3), " | ", &
            & cryst%cell(2, 1), " ", cryst%cell(2, 2), " ", cryst%cell(2, 3), " | ", &
            & cryst%cell(3, 1), " ", cryst%cell(3, 2), " ", cryst%cell(3, 3)

        do i = 1, cryst%natom
            write(102, *) cryst%name(i), cryst%xyz(i, :)
        end do
        close(unit=102)
    end subroutine

    subroutine initialize_element_map(this)
        implicit none
        class(element_map) :: this

        integer :: i

        do i = 1, 118
            this%elements(i)%number = i 
        end do

        this%elements(1)%mass = 1.008000
        this%elements(1)%symbol = "H"
        this%elements(2)%mass = 4.002602
        this%elements(2)%symbol = "He"
        this%elements(3)%mass = 6.940000
        this%elements(3)%symbol = "Li"
        this%elements(4)%mass = 9.012183
        this%elements(4)%symbol = "Be"
        this%elements(5)%mass = 10.810000
        this%elements(5)%symbol = "B"
        this%elements(6)%mass = 12.011000
        this%elements(6)%symbol = "C"
        this%elements(7)%mass = 14.007000
        this%elements(7)%symbol = "N"
        this%elements(8)%mass = 15.999000
        this%elements(8)%symbol = "O"
        this%elements(9)%mass = 18.998403
        this%elements(9)%symbol = "F"
        this%elements(10)%mass = 20.179700
        this%elements(10)%symbol = "Ne"
        this%elements(11)%mass = 22.989769
        this%elements(11)%symbol = "Na"
        this%elements(12)%mass = 24.305000
        this%elements(12)%symbol = "Mg"
        this%elements(13)%mass = 26.981538
        this%elements(13)%symbol = "Al"
        this%elements(14)%mass = 28.085000
        this%elements(14)%symbol = "Si"
        this%elements(15)%mass = 30.973762
        this%elements(15)%symbol = "P"
        this%elements(16)%mass = 32.060000
        this%elements(16)%symbol = "S"
        this%elements(17)%mass = 35.450000
        this%elements(17)%symbol = "Cl"
        this%elements(18)%mass = 39.950000
        this%elements(18)%symbol = "Ar"
        this%elements(19)%mass = 39.098300
        this%elements(19)%symbol = "K"
        this%elements(20)%mass = 40.078000
        this%elements(20)%symbol = "Ca"
        this%elements(21)%mass = 44.955908
        this%elements(21)%symbol = "Sc"
        this%elements(22)%mass = 47.867000
        this%elements(22)%symbol = "Ti"
        this%elements(23)%mass = 50.941500
        this%elements(23)%symbol = "V"
        this%elements(24)%mass = 51.996100
        this%elements(24)%symbol = "Cr"
        this%elements(25)%mass = 54.938043
        this%elements(25)%symbol = "Mn"
        this%elements(26)%mass = 55.845000
        this%elements(26)%symbol = "Fe"
        this%elements(27)%mass = 58.933154
        this%elements(27)%symbol = "Co"
        this%elements(28)%mass = 58.693400
        this%elements(28)%symbol = "Ni"
        this%elements(29)%mass = 63.546000
        this%elements(29)%symbol = "Cu"
        this%elements(30)%mass = 65.380000
        this%elements(30)%symbol = "Zn"
        this%elements(31)%mass = 69.723000
        this%elements(31)%symbol = "Ga"
        this%elements(32)%mass = 72.630000
        this%elements(32)%symbol = "Ge"
        this%elements(33)%mass = 74.921595
        this%elements(33)%symbol = "As"
        this%elements(34)%mass = 78.971000
        this%elements(34)%symbol = "Se"
        this%elements(35)%mass = 79.904000
        this%elements(35)%symbol = "Br"
        this%elements(36)%mass = 83.798000
        this%elements(36)%symbol = "Kr"
        this%elements(37)%mass = 85.467800
        this%elements(37)%symbol = "Rb"
        this%elements(38)%mass = 87.620000
        this%elements(38)%symbol = "Sr"
        this%elements(39)%mass = 88.905840
        this%elements(39)%symbol = "Y"
        this%elements(40)%mass = 91.224000
        this%elements(40)%symbol = "Zr"
        this%elements(41)%mass = 92.906370
        this%elements(41)%symbol = "Nb"
        this%elements(42)%mass = 95.950000
        this%elements(42)%symbol = "Mo"
        this%elements(43)%mass = 97.000000
        this%elements(43)%symbol = "Tc"
        this%elements(44)%mass = 101.070000
        this%elements(44)%symbol = "Ru"
        this%elements(45)%mass = 102.905490
        this%elements(45)%symbol = "Rh"
        this%elements(46)%mass = 106.420000
        this%elements(46)%symbol = "Pd"
        this%elements(47)%mass = 107.868200
        this%elements(47)%symbol = "Ag"
        this%elements(48)%mass = 112.414000
        this%elements(48)%symbol = "Cd"
        this%elements(49)%mass = 114.818000
        this%elements(49)%symbol = "In"
        this%elements(50)%mass = 118.710000
        this%elements(50)%symbol = "Sn"
        this%elements(51)%mass = 121.760000
        this%elements(51)%symbol = "Sb"
        this%elements(52)%mass = 127.600000
        this%elements(52)%symbol = "Te"
        this%elements(53)%mass = 126.904470
        this%elements(53)%symbol = "I"
        this%elements(54)%mass = 131.293000
        this%elements(54)%symbol = "Xe"
        this%elements(55)%mass = 132.905452
        this%elements(55)%symbol = "Cs"
        this%elements(56)%mass = 137.327000
        this%elements(56)%symbol = "Ba"
        this%elements(57)%mass = 138.905470
        this%elements(57)%symbol = "La"
        this%elements(58)%mass = 140.116000
        this%elements(58)%symbol = "Ce"
        this%elements(59)%mass = 140.907660
        this%elements(59)%symbol = "Pr"
        this%elements(60)%mass = 144.242000
        this%elements(60)%symbol = "Nd"
        this%elements(61)%mass = 145.000000
        this%elements(61)%symbol = "Pm"
        this%elements(62)%mass = 150.360000
        this%elements(62)%symbol = "Sm"
        this%elements(63)%mass = 151.964000
        this%elements(63)%symbol = "Eu"
        this%elements(64)%mass = 157.250000
        this%elements(64)%symbol = "Gd"
        this%elements(65)%mass = 158.925354
        this%elements(65)%symbol = "Tb"
        this%elements(66)%mass = 162.500000
        this%elements(66)%symbol = "Dy"
        this%elements(67)%mass = 164.930328
        this%elements(67)%symbol = "Ho"
        this%elements(68)%mass = 167.259000
        this%elements(68)%symbol = "Er"
        this%elements(69)%mass = 168.934218
        this%elements(69)%symbol = "Tm"
        this%elements(70)%mass = 173.045000
        this%elements(70)%symbol = "Yb"
        this%elements(71)%mass = 174.966800
        this%elements(71)%symbol = "Lu"
        this%elements(72)%mass = 178.490000
        this%elements(72)%symbol = "Hf"
        this%elements(73)%mass = 180.947880
        this%elements(73)%symbol = "Ta"
        this%elements(74)%mass = 183.840000
        this%elements(74)%symbol = "W"
        this%elements(75)%mass = 186.207000
        this%elements(75)%symbol = "Re"
        this%elements(76)%mass = 190.230000
        this%elements(76)%symbol = "Os"
        this%elements(77)%mass = 192.217000
        this%elements(77)%symbol = "Os"
        this%elements(78)%mass = 195.084000
        this%elements(78)%symbol = "Pt"
        this%elements(79)%mass = 196.966570
        this%elements(79)%symbol = "Au"
        this%elements(80)%mass = 200.592000
        this%elements(80)%symbol = "Hg"
        this%elements(81)%mass = 204.380000
        this%elements(81)%symbol = "Tl"
        this%elements(82)%mass = 207.200000
        this%elements(82)%symbol = "Pb"
        this%elements(83)%mass = 208.980400
        this%elements(83)%symbol = "Bi"
        this%elements(84)%mass = 209.000000
        this%elements(84)%symbol = "Po"
        this%elements(85)%mass = 210.000000
        this%elements(85)%symbol = "At"
        this%elements(86)%mass = 222.000000
        this%elements(86)%symbol = "Rn"
        this%elements(87)%mass = 223.000000
        this%elements(87)%symbol = "Fr"
        this%elements(88)%mass = 226.000000
        this%elements(88)%symbol = "Ra"
        this%elements(89)%mass = 227.000000
        this%elements(89)%symbol = "Ac"
        this%elements(90)%mass = 232.037700
        this%elements(90)%symbol = "Th"
        this%elements(91)%mass = 231.035880
        this%elements(91)%symbol = "Pa"
        this%elements(92)%mass = 238.028910
        this%elements(92)%symbol = "U"
        this%elements(93)%mass = 237.000000
        this%elements(93)%symbol = "Np"
        this%elements(94)%mass = 244.000000
        this%elements(94)%symbol = "Pu"
        this%elements(95)%mass = 243.000000
        this%elements(95)%symbol = "Am"
        this%elements(96)%mass = 247.000000
        this%elements(96)%symbol = "Cm"
        this%elements(97)%mass = 247.000000
        this%elements(97)%symbol = "Bk"
        this%elements(98)%mass = 251.000000
        this%elements(98)%symbol = "Cf"
        this%elements(99)%mass = 252.000000
        this%elements(99)%symbol = "Es"
        this%elements(100)%mass = 257.000000
        this%elements(100)%symbol = "Fm"
        this%elements(101)%mass = 258.000000
        this%elements(101)%symbol = "Md"
        this%elements(102)%mass = 259.000000
        this%elements(102)%symbol = "No"
        this%elements(103)%mass = 266.000000
        this%elements(103)%symbol = "Lr"
        this%elements(104)%mass = 267.000000
        this%elements(104)%symbol = "Rf"
        this%elements(105)%mass = 268.000000
        this%elements(105)%symbol = "Db"
        this%elements(106)%mass = 269.000000
        this%elements(106)%symbol = "Sg"
        this%elements(107)%mass = 270.000000
        this%elements(107)%symbol = "Bh"
        this%elements(108)%mass = 269.000000
        this%elements(108)%symbol = "Hs"
        this%elements(109)%mass = 278.000000
        this%elements(109)%symbol = "Mt"
        this%elements(110)%mass = 281.000000
        this%elements(110)%symbol = "Ds"
        this%elements(111)%mass = 282.000000
        this%elements(111)%symbol = "Rg"
        this%elements(112)%mass = 285.000000
        this%elements(112)%symbol = "Cn"
        this%elements(113)%mass = 286.000000
        this%elements(113)%symbol = "Nh"
        this%elements(114)%mass = 289.000000
        this%elements(114)%symbol = "Fl"
        this%elements(115)%mass = 290.000000
        this%elements(115)%symbol = "Mc"
        this%elements(116)%mass = 293.000000
        this%elements(116)%symbol = "Lv"
        this%elements(117)%mass = 294.000000
        this%elements(117)%symbol = "Ts"
        this%elements(118)%mass = 294.000000
        this%elements(118)%symbol = "Og"
    end subroutine

    function get_element_symbol(this, number)
        class(element_map) :: this
        character(len=2) :: get_element_symbol
        integer, intent(in) :: number
        get_element_symbol = this%elements(number)%symbol
    end function

    function get_element_number(this, symbol)
        class(element_map) :: this
        character(len=*), intent(in) :: symbol
        integer :: get_element_number
        integer :: i

        do i = 1, 118
            if (this%elements(i)%symbol == symbol) then
                get_element_number = i
                exit
            end if
        end do
    end function
end module askitf_crystal