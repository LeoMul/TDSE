module tdse_fourier
    use FFT_MOD
    use tdse
    contains 

    function get_second_derivative(psi_array,k_array) result(second_derivative)
        complex * 16 , intent(in) :: psi_array(:)
        complex * 16 :: second_derivative (size(psi_array))
        real * 8 , intent(in) :: k_array(size(psi_array))

        second_derivative = psi_array
        second_derivative = fft(second_derivative)
        second_derivative = -k_array*k_array*second_derivative
        second_derivative = ifft(second_derivative)

    end function get_second_derivative

    function get_time_derivative(psi_array,k_array,v_array) result(time_derivative)
        complex*16,intent(in) :: psi_array(:)
        real * 8,intent(in) :: k_array(:), v_array(:) 
        complex * 16 :: time_derivative(size(psi_array))

        time_derivative = get_second_derivative(psi_array,k_array)*COMPLEX(0.0_dp,-0.5_dp)
        time_derivative = time_derivative + v_array*psi_array*complex(0.0_dp,-1.0_dp)
    
    end function get_time_derivative    

    function first_time_step(psi_array,k_array,v_array,delta_t) result(newpsi)
        !https://scholars.huji.ac.il/sites/default/files/ronniekosloff/(a2At/2m)files/k19.pdf
        complex*16,intent(in) :: psi_array(:)
        real * 8,intent(in) :: k_array(:), v_array(:) ,delta_t
        complex*16 :: dpsidt (size(psi_array)),newpsi(size(psi_array)),dpsidt_2(size(psi_array))

        dpsidt = get_time_derivative(psi_array,k_array,v_array) 
        dpsidt_2  = get_time_derivative(dpsidt,k_array,v_array)
        newpsi = psi_array + dpsidt * COMPLEX(delta_t,0.0_dp) + dpsidt_2*COMPLEX(0.5_dp*delta_t*delta_t,0.0_dp)
    end function 

    function next_time_step(psi_array,psi_array_prev, k_array,v_array,delta_t) result (newpsi)
        complex*16,intent(in) :: psi_array(:),psi_array_prev(:)
        real * 8,intent(in) :: k_array(:), v_array(:) ,delta_t
        complex*16 :: dpsidt (size(psi_array)),newpsi(size(psi_array))
    
        dpsidt = get_time_derivative(psi_array,k_array,v_array)

        newpsi = psi_array_prev + dpsidt * COMPLEX(2.0_dp*delta_t,0.0_dp)


    end function next_time_step

    function propagate_wave_function_fft(initial_psi,v_array,delta_x,delta_t,num_time_steps) result (psi_matrix)
        complex*16,intent(in) :: initial_psi(:)
        real * 8,intent(in) :: delta_x, v_array(:) ,delta_t
        real * 8 :: k_array(size(v_array))
        integer :: num_time_steps , j
        complex*16 :: psi_matrix(num_time_steps,size(v_array)),next_psi(size(initial_psi))

        k_array = create_fftfreq(size(v_array),delta_x)
        psi_matrix(1,:) = initial_psi
        next_psi = first_time_step(initial_psi,k_array,v_array,delta_t)
        psi_matrix(2,:) = next_psi

        do j = 3,num_time_steps
            next_psi = next_time_step(psi_matrix(j-1,:),psi_matrix(j-2,:),k_array,v_array,delta_t)
            !print*,next_psi
            psi_matrix(j,:) = next_psi

        end do 

    end function propagate_wave_function_fft
    
    subroutine write_out_data(psi_matrix,x_array,every,delta_t,num_time_steps)
        complex*16,intent(in) :: psi_matrix(:,:) 
        real * 8, intent(in) :: x_array(:) ,delta_t
        integer,intent(in) :: every,num_time_steps
        integer :: j,jj 
        CHARACTER(LEN=20) :: FMT = "(F20.12)"
        
        open (1, file = "fourierdata.dat")
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
        write(1,*) " "

        do jj = 1,num_time_steps
            if (modulo(jj,every) == 0 ) then
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



    end subroutine 

    subroutine write_out_for_gif(psi_matrix,x_array,every,delta_t,num_time_steps)
        complex*16,intent(in) :: psi_matrix(:,:) 
        real * 8, intent(in) :: x_array(:) ,delta_t
        integer,intent(in) :: every,num_time_steps
        integer :: j,jj 
        CHARACTER(LEN=20) :: FMT = "(F20.12)"
        open (1, file = "fourierdatagif.dat")
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
    end subroutine

end module tdse_fourier