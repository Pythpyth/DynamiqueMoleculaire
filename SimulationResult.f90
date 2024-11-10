
!#####################################################real(kind=8), dimension(:), allocatable :: kinetic_energy###############################################
! SimulationResultModule
! type SimulationResult
! membres :
!      integer :: nb_time : nombre des temps de simulation
!      real(kind=8), dimension(:), allocatable :: time : tableau contenant les temps de la simulation
!      real(kind=8), dimension(:), allocatable :: kinetic_energy : tableau contenant l'�nergie cin�tique � chaque temps de la simulation
!      real(kind=8), dimension(:), allocatable :: potential_energy : tableau contenant l'�nergie potentielle � chaque temps de la simulation
!      real(kind=8), dimension(:), allocatable :: total_energy : tableau contenant l'�nergie totale � chaque temps de la simulation
!      real(kind=8), dimension(:), allocatable :: temperature : tableau contenant la temp�rature � chaque temps de la simulation
!
! subroutine :
!       terminate : d�sallocation des tableaux
!######################################################################################################


module SimulationResultModule
    use Point2DModule
    use DistanceCalculatorModule
    implicit none

    type, public :: SimulationResult
        integer :: nb_time
        real(kind=8), dimension(:), allocatable :: time
        real(kind=8), dimension(:), allocatable :: kinetic_energy
        real(kind=8), dimension(:), allocatable :: potential_energy
        real(kind=8), dimension(:), allocatable :: total_energy
        real(kind=8), dimension(:), allocatable :: temperature
    contains
        procedure, public :: terminate
    end type SimulationResult
contains
    subroutine terminate(this)
        class(SimulationResult), intent(inout) :: this
        deallocate(this%time)
        deallocate(this%kinetic_energy)
        deallocate(this%potential_energy)
        deallocate(this%total_energy)
        deallocate(this%temperature)
    end subroutine

end module SimulationResultModule



