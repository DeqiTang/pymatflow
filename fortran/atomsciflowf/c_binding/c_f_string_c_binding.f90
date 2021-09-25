module c_f_string_c_binding
    ! Usage:
    ! Reference:
    !   http://fortranwiki.org/fortran/show/c_interface_module
    use iso_c_binding, only : c_char, c_null_char, c_ptr, &
        c_associated, c_f_pointer

    implicit none

    character(len=1,kind=c_char), parameter :: NUL = c_null_char

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
        type(c_ptr), intent(in) :: c_string
        character(len=*), intent(out) :: f_string
        character(len=1,kind=C_char), dimension(:), pointer :: p_chars
        integer :: i
        if (.not. c_associated(c_string)) then
            f_string = ' '

        else
            call c_f_pointer(c_string, p_chars,[huge(0)])
            i=1
            do while(p_chars(i)/=NUL .and. i<=len(f_string))
                f_string(i:i) = p_chars(i)
                i=i+1
            end do
            if (i<len(f_string)) f_string(i:) = ' '
        end if
    end subroutine c_f_string_ptr
end module c_f_string_c_binding
