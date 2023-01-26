program main 
    use tdse_fourier
    implicit none

    real * 8, allocatable :: x_array(:),v_array(:) ,k_array(:)
    complex * 16,allocatable :: psi_array(:), psi_matrix(:,:) , psi_2(:)
    real * 8 :: delta_x , delta_t , x_last, x_first,pos,decay,strength
    integer :: N,num_time_steps,every_wave_function , i , numgauss
    !delta_x = 0.09765625_dp
    
    x_first = -20.0_dp
    x_last = 50.0_dp 
    num_time_steps = 16*32000
    N = 8192
    every_wave_function = 2000!save every 500-th wf in memory.
    allocate(psi_array(N),v_array(N),psi_matrix(num_time_steps/every_wave_function,N),k_array(N),psi_2(N),x_array(N))

    x_array = my_linspace(x_first,x_last,N)
    delta_x = x_array(2)-x_array(1)
    delta_t = delta_x*delta_x/(10.0_dp)
    print*,delta_t
    
   
    numgauss = 40
    strength = 10.0_dp
    decay = 10.0_dp
    v_array = 0.0_dp
    do i = 1,numgauss
        pos = pos + (1.0_dp / float(numgauss))  * (x_last)
        print*,"gaussian pos,",pos
        v_array = v_array + strength * exp(-decay * (x_array-pos)**2)

    end do

    psi_array = gaussian(x_array,10.0_dp,1.0_dp,0.0_dp)
    
    psi_matrix = propagate_wave_function_fft_split_operator(psi_array,v_array,delta_x,delta_t,num_time_steps,every_wave_function)
    print*,num_time_steps * delta_t
    



    call write_out_for_gif(psi_matrix,x_array,1,delta_t)
    !call write_out_data(psi_matrix,x_array,1,delta_t,num_time_steps)


end program main 