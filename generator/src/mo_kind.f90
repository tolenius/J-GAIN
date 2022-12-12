module mo_kind
    implicit none
    public
    integer, parameter ::   rk4 = kind(1.E0)
    integer, parameter ::   rk8 = kind(1.D0)
    
    integer, parameter ::   rk = rk8
    
    
    !integer(kind=2), parameter :: i2 = 1
    !integer, parameter ::   ik2 = kind(i2)
    integer(kind=4), parameter :: i4 = 1
    integer(kind=8), parameter :: i8 = 1
    integer, parameter ::   ik4 = kind(i4)
    integer, parameter ::   ik8 = kind(i8)
    
    integer, parameter ::   ik = ik4
    integer, parameter ::   dp = rk8
    
end module mo_kind