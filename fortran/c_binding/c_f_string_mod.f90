module c_f_string_mod
    ! Usage:
    ! Reference:
    !   http://fortranwiki.org/fortran/show/c_interface_module
    use iso_c_binding


    implicit none

    character(len=1,kind=C_char), parameter :: NUL = C_NULL_char

    interface c_f_string
        module procedure c_f_string_chars
        module procedure c_f_string_ptr
    end interface

    contains

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
end module c_f_string_mod
