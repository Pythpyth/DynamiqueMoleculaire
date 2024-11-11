!######################################################################################################
! GridPrinterModule
! subroutine :
! print_energy : imprime les résultats dans un fichiers
!
!######################################################################################################

module GridPrinterModule
    implicit none
contains
    subroutine print_energy(file_name, time_grid, energy_grid, grid_size)

        character(len=:), allocatable, intent(in) :: file_name
        integer, intent(in) :: grid_size
        real(kind=8), dimension (grid_size), intent(in) :: time_grid
        real(kind=8), dimension (grid_size), intent(in) :: energy_grid
        integer :: i, unit_number
        unit_number = 1

        open(unit = unit_number, file=file_name, form='formatted')

        write(unit_number, *) '#' , ' time', ' energy'
        do i =1, grid_size
            write(unit_number, *) time_grid(i) , energy_grid(i)
        end do
        close(unit_number)

    end subroutine

     subroutine print_simulation_result(file_name, time_grid, kinetic_energy, potential_energy, &
                                        total_energy, temperature, pressure, grid_size)

        character(len=:), allocatable, intent(in) :: file_name
        integer, intent(in) :: grid_size
        real(kind=8), dimension (grid_size), intent(in) :: time_grid
        real(kind=8), dimension (grid_size), intent(in) :: kinetic_energy, potential_energy, &
                                                            total_energy, temperature, pressure
        integer :: i, unit_number
        unit_number = 1

        open(unit = unit_number, file=file_name, form='formatted')

        write(unit_number, *) '#' , ' time', ' kinetic_energy', ' potential_energy', ' total_energy', &
                                ' temperature', ' pressure'
        do i =1, grid_size
            write(unit_number, *) time_grid(i) , kinetic_energy(i), potential_energy(i), &
                                 total_energy(i), temperature(i), pressure(i)
        end do

        close(unit_number)

    end subroutine

end module GridPrinterModule

