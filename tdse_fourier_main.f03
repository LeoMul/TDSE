program main 
    use tdse_fourier
    implicit none

    real * 8, allocatable :: x_array(:),v_array(:) ,k_array(:)
    complex * 16,allocatable :: psi_array(:), psi_matrix(:,:) , psi_2(:)
    real * 8 :: delta_x , delta_t , x_last, x_first
    integer :: N,num_time_steps,every_wave_function 
    !delta_x = 0.09765625_dp
    
    x_first = -20.0_dp
    x_last = 25.0_dp 
    num_time_steps = 600000
    N = 8192
    every_wave_function = 1000 !save every 500-th wf in memory.
    x_array = my_linspace(x_first,x_last,N)
    delta_x = x_array(2)-x_array(1)
    delta_t = delta_x*delta_x/(10.0_dp)
    print*,delta_t

    allocate(psi_array(N),v_array(N),psi_matrix(num_time_steps/every_wave_function,N),k_array(N),psi_2(N))
    v_array = REAL(create_square_barrier(x_array,4.0_dp,0.5_dp,50.0_dp))
    psi_array = gaussian(x_array,10.0_dp,1.0_dp,0.0_dp)
    
    psi_matrix = propagate_wave_function_fft(psi_array,v_array,delta_x,delta_t,num_time_steps,every_wave_function)
    print*,num_time_steps * delta_t
    



    call write_out_for_gif(psi_matrix,x_array,1,delta_t)
    !call write_out_data(psi_matrix,x_array,1,delta_t,num_time_steps)


end program main 