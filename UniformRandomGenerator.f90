!######################################################################################################
! UniformRandomGeneratorModule
! type UniformRandomGenerator
!       integer(kind=8) :: k : paramètre du générateur
!       integer(kind=8) :: l : paramètre du générateur
!       integer(kind=8) :: m : paramètre du générateur
!       integer(kind=8) :: x_n : valeur actuelle de xn
!
! fonction :
!       generate : le générateur prdouit des nombres aléatoires selon la règle :
!                  xn+1 = (k xn + l) mod m
!                  On en déduit des réels sur [0,1[ par :
!                  rn = xn / m
!######################################################################################################

module UniformRandomGeneratorModule
    implicit none

    type, public :: UniformRandomGenerator
        private
        integer(kind=8) :: k
        integer(kind=8) :: l
        integer(kind=8) :: m
        integer(kind=8) :: x_n
    contains
        procedure, public :: initialize, generate
    end type UniformRandomGenerator

contains
    subroutine initialize(this, k, l, m, x_0)
        class(UniformRandomGenerator), intent(inout) :: this
        integer(kind=8), intent(in) :: k, l, m, x_0
        this%k = k
        this%l = l
        this%m = m
        this%x_n = x_0
    end subroutine

    function generate(this)
        class(UniformRandomGenerator), intent(inout) :: this
        integer(kind=8) :: x_n_plus_1
        real(kind=8) :: generate
        x_n_plus_1 = modulo(this%k * this%x_n + this%l, this%m)
        this%x_n = x_n_plus_1
        generate = real(x_n_plus_1,8) / real(this%m,8)
    end function

end module UniformRandomGeneratorModule

