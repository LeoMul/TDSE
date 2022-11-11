program main 
    use tdse_fourier
    implicit none

    real * 8, allocatable :: x_array(:),v_array(:)
    complex * 16,allocatable :: psi_array(:), psi_matrix(:,:) 
    real * 8 :: delta_x , delta_t , x_last, x_first
    integer :: N,num_time_steps
    delta_x = 0.09765625_dp
    delta_t = delta_x*delta_x/(80.0_dp)
    print*,delta_t
    x_first = -25.0_dp
    x_last = 25.0_dp 
    num_time_steps = 500000


    x_array = my_arange(x_first,x_last,delta_x)
    N = size(x_array)

    allocate(psi_array(N),v_array(N),psi_matrix(num_time_steps,N))
    v_array = 0.0_dp
    psi_array = gaussian(x_array,0.0_dp,1.0_dp,0.0_dp)

    psi_matrix = propagate_wave_function_fft(psi_array,v_array,delta_x,delta_t,num_time_steps)
    !print*,psi_matrix
    call write_out_for_gif(psi_matrix,x_array,500,delta_t,num_time_steps)



end program main 