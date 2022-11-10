program main
    use tdse

    implicit none 
    CHARACTER(LEN=20) :: FMT = "(F20.12)"


    complex*16,allocatable :: m_matrix(:,:),v_array(:),identity(:,:),psi_array(:),psi_matrix(:,:),rhs(:),left_matrix(:,:),right_matrix(:,:)
    integer :: n , num_time_steps , j,jj,every
    real*8,allocatable :: x_array(:)
    real* 8 :: delta_t,delta_x ,alpha,x_0,x_N
    real* 8, parameter:: pi = 3.14159265359
    real* 8 :: strength, length, start, wavevector
    delta_t = 0.005_dp   
    x_0 = -15.0_dp  
    x_N = 15.0_dp 
    num_time_steps = 3000
    delta_x  = 0.1_dp

    alpha = 0.5_dp*delta_t/(delta_x*delta_x)
    every = 2
    x_array = my_arange(x_0,x_N,delta_x)
    n = size(x_array)
    allocate(m_matrix(n,n),v_array(n),identity(n,n),psi_array(n),psi_matrix(num_time_steps,n),rhs(n),left_matrix(n,n),right_matrix(n,n))

    identity = create_indentity_matrix(n-2) 

    start = 10.0_dp
    length = 5.0_dp 
    strength = 1.0_dp

    v_array = 0.5_dp * x_array * x_array
    m_matrix = create_m_matrix(v_array(2:n-1),COMPLEX(alpha,0.0_dp),COMPLEX(delta_t,0.0_dp))

    wavevector = 1.0_dp
    psi_array = gaussian(x_array,wavevector,1.0_dp,1.0_dp)

   
    psi_matrix(1,:) = psi_array
    print*, "solving"
    do j = 2,num_time_steps
        print*,j!,calculate_energy_integral(delta_x,v_array,psi_array)
        psi_array = schrodinger_crank_time_independent_one_step(psi_array,m_matrix,(0.0_dp,0.0_dp),(0.0_dp,0.0_dp),identity,COMPLEX(alpha,0.0_dp),delta_x)
        psi_matrix(j,:) = psi_array

    end do 
    print *, "writing out"
    open (1, file = "crank.dat")
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
        if (modulo(jj,every) == 1 ) then
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

    open (1, file = "psi_matrix_square_mod.dat")
    write(1,fmt="(1x,a)",advance = "no") "#delta_t:"

    write(1,FMT,advance = "no") delta_t
    write(1,*) " "

    do j = 1,size(x_array)
        write(1,FMT,advance = "no") x_array(j)
        write(1,fmt="(1x,a)",advance = "no") " "

    end do 
    write(1,*) " "

    do jj = 1,num_time_steps
        if (modulo(jj,every) == 1 ) then

        do j = 1,size(x_array)
            write(1,FMT,advance = "no") abs(psi_matrix(jj,j))**2 
            write(1,fmt="(1x,a)",advance = "no") " "
        end do 
        
        write(1,*) " "
         end if 
    end do 

    open (1, file = "psi_matrix_real.dat")
    write(1,fmt="(1x,a)",advance = "no") "#delta_t:"

    write(1,FMT,advance = "no") delta_t
    write(1,*) " "

    do j = 1,size(x_array)
        write(1,FMT,advance = "no") x_array(j)

    end do 
    write(1,*) " "

    do jj = 1,num_time_steps
        if (modulo(jj,every) == 1 ) then

        do j = 1,size(x_array)
            write(1,FMT,advance = "no") REALPART(psi_matrix(jj,j)) 
            write(1,fmt="(1x,a)",advance = "no") " "
        end do 
        write(1,*) " "
    end if
    end do 

    open (1, file = "psi_matrix_imag.dat")
    write(1,fmt="(1x,a)",advance = "no") "#delta_t:"

    write(1,FMT,advance = "no") delta_t
    write(1,*) " "

    do j = 1,size(x_array)
        write(1,FMT,advance = "no") x_array(j)

    end do 
    write(1,*) " "

    do jj = 1,num_time_steps
        if (modulo(jj,every) == 1 ) then

        do j = 1,size(x_array)
            write(1,FMT,advance = "no") IMAGPART(psi_matrix(jj,j)) 
            write(1,fmt="(1x,a)",advance = "no") " "
        end do 
    end if
        write(1,*) " "
    end do 

end program main