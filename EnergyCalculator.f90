!######################################################################################################
! Module EnergyCalculatorModule
! Fonctions :
!   compute_potential_energy : calcule l'énergie potentielle entre 2 particules en prenant en compte
!                              un potentiel tronqué. Si r > r_truncated, le potentiel est nul
!                              sinon on calcule le potentiel de Lennard-Jones en shiftant le potentiel de
!                              potential_shift. Le potential_shift passé en paramètre est calculé
!                              de manière à ce que le potenitel soit continu en r_truncated.
!
!   compute_sum_potential_energy : calcule l'énergie potentielle du système
!
!   compute_sum_kinetic_energy : calcule l'énergie cinétique du système
!
!
!######################################################################################################

module EnergyCalculatorModule
    use Point2DModule
    use DistanceCalculatorModule
    implicit none

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

    function compute_sum_kinetic_energy(velocities, nb_particle)
        integer(kind=4), intent(in) :: nb_particle
        type(Point2D), dimension(nb_particle), intent(in) :: velocities
        real(kind=8) :: compute_sum_kinetic_energy
        integer(kind=4) :: i

        compute_sum_kinetic_energy = 0.0d0

        do i = 1, nb_particle
             compute_sum_kinetic_energy = compute_sum_kinetic_energy + &
                                          0.5d0 * ((velocities(i)%x **2) + (velocities(i)%y ** 2))
        end do

    end function

end module EnergyCalculatorModule
