module askitf_kinds
    
    ! fortran 2008
    !use, intrinsic :: iso_fortran_env
    
    implicit none

    ! fortran 2008
    !integer, parameter :: sp = REAL32
    !integer, parameter :: dp = REAL64
    !integer, parameter :: qp = REAL128
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: qp = selected_real_kind(33, 4931)
    
end module askitf_kinds