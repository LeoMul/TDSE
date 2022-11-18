program main 
    use tdse_fourier_2d
    use tdse_fourier
    use tdse
    implicit none
    integer,parameter :: Nx = 256, Ny = 256
    real * 8 :: x_array_temp(Nx),y_array_temp(Ny), k_x(Nx),k_y(Ny)
    real * 8 :: delta_x, delta_y,first_x,final_x ,delta_t,first_y,final_y,wavenumber
    complex * 16,allocatable :: psi_array(:,:),final_psi(:,:,:)
    real * 8,allocatable :: k_matrix(:,:),v_array(:,:)
    integer :: num_time_steps ,every 
    first_x = -5_dp 
    final_x = 5_dp 
    first_y = -5_dp
    final_y = 5_dp
    x_array_temp = my_linspace(first_x,final_x,Nx)
    y_array_temp = my_linspace(first_y,final_y,Ny) 

    delta_x = x_array_temp(2) - x_array_temp(1)
    delta_y = y_array_temp(2) - y_array_temp(1)

    delta_t = (delta_x**2 + delta_y**2)/50.0_dp

    k_x =  6.28318530718_dp*create_fftfreq(Nx,delta_x)
    k_y =  6.28318530718_dp*create_fftfreq(Ny,delta_y)

    wavenumber = 5.0_dp
    psi_array = gaussian_2d(x_array_temp,y_array_temp,1.0_dp,wavenumber)
    k_matrix = create_k_matrix(k_x,k_y)

    v_array = 0.0_dp * k_matrix

    !v_array = create_circular_barrier_potential(x_array_temp,y_array_temp,10000.0_dp,3.0_dp,0.0_dp,1.0_dp)
    !print*,v_array

    !psi_array = get_laplacian(psi_array,k_matrix)

    print *,size(psi_array)
    print*,size(psi_array,1)
    print*,size(psi_array,2)
    num_time_steps = 10000
    every = 50

    !print*,k_x
    !print*,k_y
    !print*,k_matrix
    final_psi = propagate_wave_function_fft_2d(psi_array,v_array,delta_t,k_matrix,num_time_steps,every)
    

    !call write_out_dim2array(psi_array,x_array_temp,y_array_temp,"initipsi.dat")

    call write_out_3index_for_gif(final_psi,x_array_temp,y_array_temp,every,delta_t)

end program 