!######################################################################################################
! ObservableCalculatorModule
! Fonctions :
!   compute_mean_observable_from_t_ref : calcul la moyenne temporelle d'un observable à partir du temps t_ref
!
!   compute_stdev_observable_from_t_ref : calcul de l'écart type temporelle d'un observable à partir du temps t_ref
!
!  compute_number_observable_from_t_ref : calcule le nombre d'observale depuis t_ref
!
!######################################################################################################


module ObservableCalculatorModule
    implicit none


contains
    function compute_mean_observable_from_t_ref(t_ref, times, observable, nb_times)
        integer, intent(in) :: nb_times
        real(kind=8), intent(in) :: t_ref
        real(kind=8), dimension(nb_times), intent(in) :: times, observable
        real(kind=8) :: compute_mean_observable_from_t_ref

        integer :: i, nb_times_for_mean

        nb_times_for_mean = 0
        compute_mean_observable_from_t_ref = 0.0d0

        do i = 1, nb_times

            if (times(i).ge.t_ref) then
                nb_times_for_mean = nb_times_for_mean + 1
                compute_mean_observable_from_t_ref = compute_mean_observable_from_t_ref + observable(i)
            end if

        end do

       compute_mean_observable_from_t_ref  = compute_mean_observable_from_t_ref / real(nb_times_for_mean, 8)
    end function

    function compute_stdev_observable_from_t_ref(t_ref, times, observable, nb_times)
        integer, intent(in) :: nb_times
        real(kind=8), intent(in) :: t_ref
        real(kind=8), dimension(nb_times), intent(in) :: times, observable
        real(kind=8) :: compute_stdev_observable_from_t_ref
        real(kind=8) :: mean, mean_square, std_dev

        integer :: i, nb_times_for_mean

        nb_times_for_mean = 0
        mean = 0.0d0
        mean_square = 0.0d0

        do i = 1, nb_times

            if (times(i).ge.t_ref) then
                nb_times_for_mean = nb_times_for_mean + 1
                mean = mean + observable(i)
                mean_square = mean_square + observable(i)**2
            end if

        end do

       std_dev  =  sqrt( mean_square / real(nb_times_for_mean, 8) -  ((mean / real(nb_times_for_mean, 8))**2) )
      compute_stdev_observable_from_t_ref = std_dev
    end function

    function compute_number_observable_from_t_ref(t_ref, times, nb_times)
        integer, intent(in) :: nb_times
        real(kind=8), intent(in) :: t_ref
        real(kind=8), dimension(nb_times), intent(in) :: times
        integer :: compute_number_observable_from_t_ref

        integer :: i

        compute_number_observable_from_t_ref = 0

        do i = 1, nb_times

            if (times(i).ge.t_ref) then
                compute_number_observable_from_t_ref = compute_number_observable_from_t_ref + 1
            end if

        end do

    end function


end module ObservableCalculatorModule

