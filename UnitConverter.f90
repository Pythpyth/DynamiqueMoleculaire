
!######################################################################################################
! UnitConverterModule
! Fonctions :
!   convert_reduced_unit_to_energy_unit : E(unité énergie) = E(unité réduite) * epsilon Lennard-Jones(unité énergie)
!
!   convert_reduced_unit_to_kelvin : T(K) = T(unité réduite) * epsilon Lennard-Jones(unité atomique) / kb(unité atomique par K)
!
!   convert_reduced_unit_to_second : t(unité atomic) = t(unité réduite) * racine(masse particule(unité atomique) *
!                                                                       sigma Lennard-Jones(unité atomique)^2 /
!                                                                       epsilon Lennard-Jones(unité atomique)^2)
!                                    t(fs) =   t(unité atomic)  * atomic_time_to_fs
!
!
!
!######################################################################################################

module UnitConverterModule
    use Point2DModule
    use DistanceCalculatorModule
    implicit none

contains

    function convert_reduced_unit_to_energy_unit(energy_list, list_size, epsilon_LJ)
        integer, intent(in) :: list_size
        real(kind=8), dimension (list_size), intent(in) :: energy_list
        real(kind=8), dimension (list_size) :: convert_reduced_unit_to_energy_unit
        real(kind=8) :: epsilon_LJ

        integer :: i
        do i = 1 , list_size
            convert_reduced_unit_to_energy_unit(i) = energy_list(i) * epsilon_LJ
        end do
    end function

    function convert_list_reduced_unit_to_kelvin(temperature_list, list_size, epsilon_LJ_atomic_unit, &
                                                 kb_atomic_by_K)
        integer, intent(in) :: list_size
        real(kind=8), dimension (list_size), intent(in) :: temperature_list
        real(kind=8), dimension (list_size) :: convert_list_reduced_unit_to_kelvin
        real(kind=8) :: epsilon_LJ_atomic_unit, kb_atomic_by_K
        integer :: i

        do i = 1 , list_size

            convert_list_reduced_unit_to_kelvin(i) = convert_reduced_unit_to_kelvin(temperature_list(i), &
                                                                                epsilon_LJ_atomic_unit , kb_atomic_by_K)
        end do

    end function

    function convert_reduced_unit_to_kelvin(temperature, epsilon_LJ_atomic_unit, kb_atomic_by_K)
        real(kind=8) convert_reduced_unit_to_kelvin
        real(kind=8) :: temperature, epsilon_LJ_atomic_unit, kb_atomic_by_K

        convert_reduced_unit_to_kelvin = temperature * epsilon_LJ_atomic_unit / kb_atomic_by_K

    end function


    function convert_list_reduced_unit_to_second(time_list, list_size, particle_mass_atomic_unit, &
                                                 epsilon_LJ_atomic_unit, sigma_LJ_atomic_unit)
        integer, intent(in) :: list_size
        real(kind=8), dimension (list_size), intent(in) :: time_list
        real(kind=8), dimension (list_size) :: convert_list_reduced_unit_to_second
        real(kind=8) :: particle_mass_atomic_unit, epsilon_LJ_atomic_unit, sigma_LJ_atomic_unit
        integer :: i

        do i = 1 , list_size

            convert_list_reduced_unit_to_second(i) = &
            convert_reduced_unit_to_second(time_list(i), particle_mass_atomic_unit, &
                                           epsilon_LJ_atomic_unit, sigma_LJ_atomic_unit)
        end do

    end function

    function convert_reduced_unit_to_second(time, particle_mass_atomic_unit, epsilon_LJ_atomic_unit, &
                                            sigma_LJ_atomic_unit)
        real(kind=8) :: convert_reduced_unit_to_second
        real(kind=8) :: time, particle_mass_atomic_unit, epsilon_LJ_atomic_unit, sigma_LJ_atomic_unit
        real(kind=8), parameter :: atomic_time_to_fs = 0.024188843265864


        convert_reduced_unit_to_second = time * sqrt(particle_mass_atomic_unit* (sigma_LJ_atomic_unit**2) /&
                                                        epsilon_LJ_atomic_unit) * atomic_time_to_fs

    end function



end module UnitConverterModule
