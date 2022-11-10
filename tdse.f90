program main
    !gfortran -fbacktrace -Wall -fcheck=all tdse_mod.f90 tdse.f90 -ffree-line-length-none -o int_check
    use tdse
    
    implicit none
    CHARACTER(LEN=20) :: FMT = "(F20.12)"

    complex*16,allocatable :: psi_array(:),psi_matrix(:,:)
    real*8 :: alpha, delta_x , delta_t , x_0, x_N ,t_0, t_final ,k
    real*8,allocatable :: x_array(:),V_array(:)
    integer :: j,jj,num_time_steps

    x_0 = -10.0_dp
    x_N = -x_0

    num_time_steps = 200000
    delta_x = 0.05_dp
    delta_t = delta_x**2/100.0_dp
    print*,"delta_t",delta_t
    alpha = delta_t/delta_x**2
    t_0 = 0.0_dp 
    t_final = num_time_steps*delta_t

    k = 5.0_dp
    
    x_array = my_arange(x_0,X_N,delta_x)
    allocate(psi_array(size(x_array)))

    psi_array = gaussian(x_array,k,1.0_dp,0.0_dp)
    psi_array(1) = (0.0_dp,0.0_dp)
    psi_array(size(psi_array))= (0.0_dp,0.0_dp)
    V_array = create_v_array(x_array)
    !print*,V_array
    allocate(psi_matrix(num_time_steps,size(x_array)))

    psi_matrix(1,:) = psi_array

    do j = 2,num_time_steps
        psi_array = calculate_psi_next_time_step(psi_array,(0.0_dp,0.0_dp),(0.0_dp,0.0_dp),V_array,alpha,delta_t,delta_x)
        psi_matrix(j,:) = psi_array 
        print*,j   
    end do 
    !print*,psi_matrix

    open (1, file = "file_name_1.dat")
    J = 1 !every jth point
    write(1,fmt="(1x,a)",advance = "no") "    "


    write(1,fmt="(1x,a)",advance = "no") "#t"
    write(1,fmt="(1x,a)",advance = "no") "      "
    write(1,fmt="(1x,a)",advance = "no") "#x"
    write(1,fmt="(1x,a)",advance = "no") "      "
    write(1,fmt="(1x,a)",advance = "no") "#Re(psi)"
    write(1,fmt="(1x,a)",advance = "no") "  "
    write(1,fmt="(1x,a)",advance = "no") "#Im(psi)"
    write(1,fmt="(1x,a)",advance = "no") "  "
    write(1,fmt="(1x,a)",advance = "no") "#abs(psi)^2"
    write(1,fmt="(1x,a)",advance = "no") "  234"
    write(1,*) " "
    
    do jj = 1,num_time_steps
        if (modulo(jj,1000) == 1 ) then
    do j = 1,size(x_array)
                write(1,FMT,advance = "no") (jj-1)*delta_t
                write(1,fmt="(1x,a)",advance = "no") " "
                write(1,FMT,advance = "no") x_array(j)
                write(1,fmt="(1x,a)",advance = "no") " "
                write(1,FMT,advance = "no") REALPART(psi_matrix(jj,j)) 
                write(1,fmt="(1x,a)",advance = "no") " "
                write(1,FMT,advance = "no") IMAGPART(psi_matrix(jj,j)) 
                write(1,fmt="(1x,a)",advance = "no") " "
                write(1,FMT,advance = "no") abs(psi_matrix(jj,j))**2 
                write(1,*) " "

    end do 
end if
    end do 

end program 



