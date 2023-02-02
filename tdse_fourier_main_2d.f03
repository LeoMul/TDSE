program main 
    use tdse_fourier_2d
    use tdse_fourier
    use tdse
    implicit none
    integer,parameter :: Nx = 512, Ny = 512
    real * 8 :: x_array_temp(Nx),y_array_temp(Ny), k_x(Nx),k_y(Ny)
    real * 8 :: delta_x, delta_y,first_x,final_x ,delta_t,first_y,final_y,wavenumber

    real * 8 :: strength, width, spacing , packetwidth ,xthresh

    complex * 16,allocatable :: psi_array(:,:),final_psi(:,:,:)
    real * 8,allocatable :: k_matrix(:,:),v_array(:,:)
    integer :: num_time_steps ,every 

    first_x = -15.0_dp 
    final_x = 20.0_dp 
    first_y = -2.5_dp
    final_y = 2.5_dp

    x_array_temp = my_linspace(first_x,final_x,Nx)
    y_array_temp = my_linspace(first_y,final_y,Ny) 

    delta_x = x_array_temp(2) - x_array_temp(1)
    delta_y = y_array_temp(2) - y_array_temp(1)

    delta_t = (delta_x**2 + delta_y**2)/400.0_dp

    num_time_steps = 10000
    every = 400

    print*,"delta x ",delta_x
    print*,"delta y ",delta_y
    print*,"delta t ",delta_t
    print*,"total time",delta_t*num_time_steps

    k_x =  6.28318530718_dp*create_fftfreq(Nx,delta_x)
    k_y =  6.28318530718_dp*create_fftfreq(Ny,delta_y)
    allocate(v_array(Nx,Ny))

    wavenumber = 10.0_dp
    packetwidth = 2.0_dp
    psi_array = gaussian_2d(x_array_temp,y_array_temp,packetwidth,wavenumber)
    k_matrix = create_k_matrix(k_x,k_y)
    v_array = 0.0_dp

    spacing = 2.5_dp 
    strength = 100.0_dp
    width = 0.25_dp
    xthresh = 1.0_dp

    v_array = create_gaussian_crystal_potential(x_array_temp,y_array_temp,xthresh,spacing,strength,width)

    call write_out_potential(v_array,strength,spacing,width)

    print*,"Potential Creation complete"
    !v_array = create_circular_barrier_potential(x_array_temp,y_array_temp,10000.0_dp,3.0_dp,0.0_dp,1.0_dp)
    !print*,v_array

    !psi_array = get_laplacian(psi_array,k_matrix)

    print *,"total elements per time step",size(psi_array)
    print*,"Nx",size(psi_array,1)
    print*,"Ny",size(psi_array,2)
    print*,"max kx", 6.28318530718_dp * real(Nx/2 - 1) / (Nx*delta_x)
    print*,"max ky",6.28318530718_dp * real(Ny/2 - 1) / (Ny*delta_x)


    !print*,k_x
    !print*,k_y
    !print*,k_matrix
    final_psi = propagate_wave_function_fft_2d_split_operator(psi_array,v_array,k_matrix,delta_t,num_time_steps,every)
    

    !call write_out_dim2array(psi_array,x_array_temp,y_array_temp,"initipsi.dat")

    call write_out_3index_for_gif_gaussian(final_psi,x_array_temp,y_array_temp,every,delta_t,width,strength,spacing,packetwidth,wavenumber)

end program 