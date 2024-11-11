
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


    function convert_list_reduced_unit_to_fs(time_list, list_size, particle_mass_atomic_unit, &
                                                 epsilon_LJ_atomic_unit, sigma_LJ_atomic_unit)
        integer, intent(in) :: list_size
        real(kind=8), dimension (list_size), intent(in) :: time_list
        real(kind=8), dimension (list_size) :: convert_list_reduced_unit_to_fs
        real(kind=8) :: particle_mass_atomic_unit, epsilon_LJ_atomic_unit, sigma_LJ_atomic_unit
        integer :: i

        do i = 1 , list_size

            convert_list_reduced_unit_to_fs(i) = &
            convert_reduced_unit_to_fs(time_list(i), particle_mass_atomic_unit, &
                                           epsilon_LJ_atomic_unit, sigma_LJ_atomic_unit)
        end do

    end function

    function convert_reduced_unit_to_fs(time, particle_mass_atomic_unit, epsilon_LJ_atomic_unit, &
                                            sigma_LJ_atomic_unit)
        real(kind=8) :: convert_reduced_unit_to_fs
        real(kind=8) :: time, particle_mass_atomic_unit, epsilon_LJ_atomic_unit, sigma_LJ_atomic_unit
        real(kind=8), parameter :: atomic_time_to_fs = 0.024188843265864d0


        convert_reduced_unit_to_fs = time * sqrt(particle_mass_atomic_unit* (sigma_LJ_atomic_unit**2) /&
                                                        epsilon_LJ_atomic_unit) * atomic_time_to_fs

    end function

     function convert_list_reduced_unit_to_newton_by_meter(pressure, list_size, &
                                                            epsilon_LJ_atomic_unit, sigma_LJ_atomic_unit)
        integer, intent(in) :: list_size
        real(kind=8), dimension(list_size), intent(in) :: pressure
        real(kind=8), intent(in) :: epsilon_LJ_atomic_unit, sigma_LJ_atomic_unit
        real(kind=8), dimension(list_size) :: convert_list_reduced_unit_to_newton_by_meter
        integer :: i

        do i = 1 , list_size

            convert_list_reduced_unit_to_newton_by_meter(i) = &
            convert_reduced_unit_to_newton_by_meter(pressure(i), epsilon_LJ_atomic_unit, sigma_LJ_atomic_unit)
        end do
    end function

    function convert_reduced_unit_to_newton_by_meter(pressure, epsilon_LJ_atomic_unit, sigma_LJ_atomic_unit)
        real(kind=8), intent(in) :: pressure, epsilon_LJ_atomic_unit, sigma_LJ_atomic_unit
        real(kind=8) :: convert_reduced_unit_to_newton_by_meter
        !newton_by_meter = joule_by_meter_square
        real(kind=8), parameter :: hartree_by_bohr_square_to_newton_by_meter = 1556.89310492614d0


        convert_reduced_unit_to_newton_by_meter = hartree_by_bohr_square_to_newton_by_meter * &
                                                pressure * epsilon_LJ_atomic_unit / (sigma_LJ_atomic_unit**2)

    end function



end module UnitConverterModule
