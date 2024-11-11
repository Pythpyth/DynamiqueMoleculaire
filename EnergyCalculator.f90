!######################################################################################################
! Module EnergyCalculatorModule
! Fonctions :
!   compute_potential_energy : calcul l'énergie potentielle entre 2 particules en prenant en compte
!                              un potentiel tronqué. Si r > r_truncated, le potentiel est nul
!                              sinon on calcule le potentiel de Lennard-Jones en shiftant le potentiel de
!                              potential_shift. Le potential_shift passé en paramètre est calculé
!                              de manière à ce que le potenitel soit continu en r_truncated.
!
!   compute_sum_potential_energy : calcul l'énergie potentielle du système
!
!   compute_sum_kinetic_energy : calcul l'énergie cinétique du système
!
!   compute_temperature : calcul de la température Ec = N * kb * T en 2D en unité réduite kb = 1
!                                                  T = Ec /N
!
!
!   compute_pressure : calcul de la pression : P = N * kb * T /Aire + 1/(2*Aire) * somme(force . distance) en 2D en unité réduite kb = 1
!                                              P = N * T /Aire + 1/(2*Aire) * somme(force . distance)
!   compute_potential_energy_tail_correction : calcul de la correction en 2 d de queue de l'énergie potentielle
!                                              r_truncated et box_size sont passés sans dimension
!                                              on a: r_truncated_dimensioné = r_truncated * sigma_LJ
!                                                    box_size_dimensioné = box_size * sigma_LG
!
!
!   compute_potential_energy_shift_correction : calcul de la correction en 2 d liée au shift de l'énergie potentielle
!                                               correction = 0.5 * nb_particle * <nc> * potential_shift
!                                               ou <nc> représente le nombre moyen de particules dans le volume
!                                               correspondant au rayon de coupure rc
!                                               on estime : <rc> = aire cercle de rayon rc * densité de particule dans la boite
!
!
!   compute_pressure_tail_correction : calcul de la correction en 2 d de queue de la pression
!                                      r_truncated et box_size sont passés sans dimension
!                                      on a: r_truncated_dimensioné = r_truncated * sigma_LJ
!                                      box_size_dimensioné = box_size * sigma_LJ
!
!######################################################################################################

module EnergyCalculatorModule
    use Point2DModule
    use DistanceCalculatorModule
    use UnitConverterModule
    use ForceCalculatorModule, only : ForceCalculator
    implicit none
    real(kind=8), parameter:: pi = 3.1415926535897932d0

contains

    function compute_potential_energy(position_reference, position_2, &
                                      box_size, half_box_size, r_truncated, potential_shift)
        type(Point2D), intent(in):: position_reference, position_2
        real(kind=8), intent(in):: box_size, half_box_size, r_truncated, potential_shift
        type(Point2D) :: closest_image_position_2
        real(kind=8) :: compute_potential_energy
        real(kind=8) :: distance

        closest_image_position_2 = compute_closest_image_position_2(&
        position_reference, position_2, half_box_size, box_size)

        distance = compute_distance(position_reference, closest_image_position_2)

        if (distance.gt.(r_truncated)) then
            compute_potential_energy = 0.0d0
        else
            compute_potential_energy = 4.0d0 * ( distance ** (-12) - distance ** (-6)) - potential_shift
        end if

    end function

    function compute_sum_potential_energy(positions, nb_particle, box_size, half_box_size, &
                                          r_truncated, potential_shift)
        integer(kind=4), intent(in) :: nb_particle
        real(kind=8), intent(in):: box_size, half_box_size, r_truncated, potential_shift
        type(Point2D), dimension(nb_particle), intent(in) :: positions
        real(kind=8) :: compute_sum_potential_energy
        integer(kind=4) :: i,j

        compute_sum_potential_energy = 0.0d0

        do i = 1, nb_particle

             do j = i + 1, nb_particle
                 compute_sum_potential_energy = compute_sum_potential_energy &
                 + compute_potential_energy(positions(i), positions(j), box_size, half_box_size, &
                                            r_truncated, potential_shift)
             end do

        end do

    end function


    function compute_potential_energy_tail_correction(nb_particle, box_size, r_truncated, &
                                                      epsilon_LJ)

        integer(kind=4), intent(in) :: nb_particle
        real(kind=8), intent(in):: box_size, r_truncated, epsilon_LJ

        real(kind=8) :: compute_potential_energy_tail_correction
        real(kind=8) :: density

        density = real(nb_particle, 8) /(box_size * box_size)

        compute_potential_energy_tail_correction = 4.0d0 * pi * real(nb_particle, 8) * density * epsilon_LJ * &
                                                   ( 0.1d0 * (r_truncated **(-10)) - 0.25d0 * (r_truncated **(-4)) )

    end function

    function compute_potential_energy_shift_correction(nb_particle, box_size, r_truncated, &
                                                       potential_shift, epsilon_LJ)

        integer(kind=4), intent(in) :: nb_particle
        real(kind=8), intent(in):: box_size, r_truncated, potential_shift, epsilon_LJ

        real(kind=8) :: compute_potential_energy_shift_correction
        real(kind=8) :: density

        density = real(nb_particle, 8) /(box_size * box_size)

        compute_potential_energy_shift_correction = 0.5d0 * real(nb_particle,8) * density * pi * (r_truncated**2) &
                                                    * potential_shift * epsilon_LJ
    end function


    function compute_sum_kinetic_energy(velocities, nb_particle)
        integer(kind=4), intent(in) :: nb_particle
        type(Point2D), dimension(nb_particle), intent(in) :: velocities
        real(kind=8) :: compute_sum_kinetic_energy
        integer(kind=4) :: i

        compute_sum_kinetic_energy = 0.0d0

        do i = 1, nb_particle
             compute_sum_kinetic_energy = compute_sum_kinetic_energy + 0.5d0 * ((velocities(i)%x **2) + (velocities(i)%y ** 2))
        end do

        compute_sum_kinetic_energy = compute_sum_kinetic_energy

    end function

    function compute_temperature(kinetic_energy, nb_particle)
        integer(kind=4), intent(in) :: nb_particle
        real(kind=8), intent(in) :: kinetic_energy
        real(kind=8) :: compute_temperature

        compute_temperature = kinetic_energy / (real(nb_particle,8))

    end function

    function compute_pressure(positions, nb_particle, box_size, temperature, force_calculator)
        integer(kind=4), intent(in) :: nb_particle
        type(ForceCalculator), intent(in)  :: force_calculator
        real(kind=8), intent(in):: box_size, temperature
        type(Point2D), dimension(nb_particle), intent(in) :: positions
        real(kind=8) :: entropy_component, area
        real(kind=8) :: compute_pressure
        integer(kind=4) :: i,j

        area= box_size * box_size


        compute_pressure = 0.0d0
        do i = 1, nb_particle

             do j = i + 1, nb_particle
                 compute_pressure = compute_pressure + &
                 force_calculator%compute_force_scalar_distance_on_reference_particle(positions(i), positions(j))
             end do

        end do

        compute_pressure = compute_pressure * 0.5d0 / area

        entropy_component = real(nb_particle,8) * temperature / area

        compute_pressure = compute_pressure + entropy_component

    end function

    function compute_pressure_tail_correction(nb_particle, box_size, r_truncated, epsilon_LJ_atomic_unit, sigma_LJ_atomic_unit)

        integer(kind=4), intent(in) :: nb_particle
        real(kind=8), intent(in):: box_size, r_truncated, epsilon_LJ_atomic_unit, sigma_LJ_atomic_unit
        real(kind=8) :: compute_pressure_tail_correction
        real(kind=8) :: density

        density = real(nb_particle, 8) /(box_size * box_size)

        compute_pressure_tail_correction = 12.0d0 * pi * (density**2) * &
                                            (0.2d0 * (r_truncated**(-10)) - 0.25d0 * (r_truncated**(-4)))

        compute_pressure_tail_correction = convert_reduced_unit_to_newton_by_meter(compute_pressure_tail_correction , &
                                                                       epsilon_LJ_atomic_unit, sigma_LJ_atomic_unit)


    end function

end module EnergyCalculatorModule
