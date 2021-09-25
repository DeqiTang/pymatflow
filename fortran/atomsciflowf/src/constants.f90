module asflowf_constants
    use asflowf_kinds, only : dp

    implicit none
    real(kind=dp), parameter :: bohr_to_angstrom = 0.529177249
    real(kind=dp), parameter :: pi = 3.141592653
    real(kind=dp), parameter :: ry_to_ev = 13.6056923
    real(kind=dp), parameter :: ha_to_ev = 27.211324570273
    
end module asflowf_constants