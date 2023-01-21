module tdse_fourier_2d
    use FFT_MOD
    use tdse
    contains

    function get_laplacian(psi_array,k_matrix) result(second_derivative)
        complex * 16 , intent(in) :: psi_array(:,:)
        complex * 16 :: second_derivative (size(psi_array,1),size(psi_array,2)) !check this
        real * 8,intent(in) :: k_matrix(:,:)


        second_derivative = psi_array
        second_derivative = fft_2d(second_derivative)

        second_derivative = -k_matrix*second_derivative
        second_derivative = ifft_2d(second_derivative)

    end function get_laplacian

    function create_k_matrix (k_x,k_y) result (k_squared_matrix)
        real *8 ,intent(in):: k_x(:),k_y(:) 
        integer :: n_x,n_y ,i,j
        real * 8:: k_x_squared(size(k_x)),k_y_squared(size(k_y))
        real * 8 :: k_squared_matrix(size(k_x),size(k_y))
        n_x = size(k_x)
        n_y = size(k_y)
        k_x_squared = k_x * k_x
        k_y_squared = k_y * k_y 
        !this may not be correct yet for arbitrary size. 
        do i = 1,n_x
            do j = 1,n_y 
                k_squared_matrix(i,j) = k_x_squared(i) + k_y_squared(j)
            end do 
        end do 
    end function

    function get_time_derivative(psi_array,k_array,v_array) result(time_derivative)
        complex*16,intent(in) :: psi_array(:,:)
        real * 8,intent(in) :: k_array(:,:), v_array(:,:) 
        complex * 16 :: time_derivative(size(psi_array,1),size(psi_array,2))

        time_derivative = get_laplacian(psi_array,k_array)*COMPLEX(0.0_dp,0.5_dp)
        time_derivative = time_derivative + v_array*psi_array*complex(0.0_dp,-1.0_dp)
    
    end function get_time_derivative    

    function first_time_step(psi_array,k_array,v_array,delta_t) result(newpsi)
        !https://scholars.huji.ac.il/sites/default/files/ronniekosloff/(a2At/2m)files/k19.pdf
        complex*16,intent(in) :: psi_array(:,:)
        real * 8,intent(in) :: k_array(:,:), v_array(:,:) ,delta_t
        complex*16 :: dpsidt (size(psi_array,1),size(psi_array,2)),newpsi(size(psi_array,1),size(psi_array,2)),dpsidt_2(size(psi_array,1),size(psi_array,2))

        dpsidt = get_time_derivative(psi_array,k_array,v_array) 
        dpsidt_2  = get_time_derivative(dpsidt,k_array,v_array)
        newpsi = psi_array + dpsidt * COMPLEX(delta_t,0.0_dp) + dpsidt_2*COMPLEX(0.5_dp*delta_t*delta_t,0.0_dp)
    end function 

    function next_time_step(psi_array,psi_array_prev, k_array,v_array,delta_t) result (newpsi)
        complex*16,intent(in) :: psi_array(:,:),psi_array_prev(:,:)
        real * 8,intent(in) :: k_array(:,:), v_array(:,:) ,delta_t
        complex*16 :: dpsidt (size(psi_array,1),size(psi_array,2)),newpsi(size(psi_array,1),size(psi_array,2))
    
        dpsidt = get_time_derivative(psi_array,k_array,v_array)

        newpsi = psi_array_prev + dpsidt * COMPLEX(2.0_dp*delta_t,0.0_dp)


    end function next_time_step


    function propagate_wave_function_fft_2d(initial_psi,v_array,delta_t,k_matrix,num_time_steps,every) result (psi_3index)
        complex*16,intent(in) :: initial_psi(:,:)
        
        real *8,intent(in) :: k_matrix(:,:),delta_t ,v_array(:,:)
        integer, intent(in) ::num_time_steps ,every 
        complex*16 :: next_psi(size(initial_psi,1),size(initial_psi,2)) ,prev_psi(size(initial_psi,1),size(initial_psi,2)) ,current_psi(size(initial_psi,1),size(initial_psi,2))
        integer :: j 
        complex * 16 :: psi_3index(num_time_steps/every+1,size(initial_psi,1),size(initial_psi,2))
        prev_psi = initial_psi
        current_psi = first_time_step(initial_psi,k_matrix,v_array,delta_t)
        psi_3index(1,:,:) = current_psi
        do j = 1,num_time_steps
            next_psi = next_time_step(current_psi,prev_psi,k_matrix,v_array,delta_t)
            prev_psi = current_psi
            current_psi = next_psi
            !print*,shape(current_psi)
            if ((modulo(j,every) == 0) .or. (every == 1)) then
                print*,j
 
                psi_3index(j/every+1,:,:) = next_psi
            end if 
            

        end do
    end function 

    function propagate_wave_function_fft_2d_split_operator(initial_psi,v_array,k_matrix,delta_t,num_time_steps,every) result (psi_3index)

        complex*16,intent(in) :: initial_psi(:,:)
        integer ,intent(in) :: every
        real * 8,intent(in) :: v_array(:,:) ,delta_t,k_matrix(:,:)
        integer :: num_time_steps , j
        complex * 16::psi(size(initial_psi,1),size(initial_psi,2)),exptarray(size(v_array,1),size(v_array,2)),expv_array(size(v_array,1),size(v_array,2))
        complex * 16 :: psi_3index(num_time_steps/every+1,size(initial_psi,1),size(initial_psi,2))

        psi_3index(1,:,:) = initial_psi

        exptarray = exp(complex(0.0_dp,-0.5_dp*delta_t)*k_matrix*k_matrix)
        expv_array= exp(complex(0.0_dp,-0.5_dp*delta_t)*v_array)
        psi = initial_psi

        do j = 2,num_time_steps
            psi = expv_array*psi
            psi = fft_2d(psi)
            psi = exptarray * psi 
            psi = expv_array * ifft_2d(psi)
            if ((modulo(j,every) == 0) .or. (every == 1)) then
                print*,j
 
                psi_3index(j/every+1,:,:) = psi
            end if 
        end do 

    end function propagate_wave_function_fft_2d_split_operator


    function create_circular_barrier_potential(x_array,y_array,strength,x_pos,y_pos,radius) result(v_array)
        real * 8,intent(in) :: x_array(:),y_array(:)
        real * 8,intent(in) :: strength, x_pos, y_pos, radius
        real *8 :: v_array(size(x_array),size(y_array)) ,dist,radsq
        integer :: i ,j 
        radsq = radius * radius
        do i = 1,size(x_array)
            do j = 1,size(y_array)
                dist = (x_array(i)-x_pos)**2 + (y_array(j)-y_pos)**2 
                if (dist .LE. radsq) then 
                    v_array(i,j) = strength
                else 
                    !print*,20
                    v_array(i,j) = 0.0_dp
                end if 
            end do 
        end do
    end function 

    function create_gaussian_crystal_potential(x_array,y_array,xthresh,ystart,spacing,strength,width) result (v_array)
        real * 8,intent(in) :: x_array(:),y_array(:)
        real * 8,intent(in) :: strength, xthresh,spacing,ystart,width
        real *8 :: v_array(size(x_array),size(y_array)),currentx,currenty
        integer :: i,j,Nx,Ny

        v_array = 0.0_dp

        Ny = NINT((y_array(size(y_array)) - y_array(1)) /spacing )
        Nx = NINT((x_array(size(x_array)) - xthresh) /spacing ) - 1

        print*,"Requested" ,Nx,Ny,"gaussians"

        currentx = xthresh + spacing

        do i = 1,Nx
            currenty = ystart
            do j = 1,Ny 
                v_array = v_array + strength*gaussian_2d_potential(x_array,y_array,width,currentx,currenty)
                currenty = currenty + spacing
            end do 
            currentx = currentx + spacing
        end do 

    end function


    subroutine write_out_dim2array(psi_matrix,x_array,y_array,name)
        complex*16,intent(in) :: psi_matrix(:,:) 
        real * 8, intent(in) :: x_array(:) ,y_array(:)
        character(LEN = 12),intent(in) :: name
        integer :: j,jj 
        CHARACTER(LEN=20) :: FMT = "(F20.12)"
        
        open (1, file = name)
        write(1,fmt="(1x,a)",advance = "no") "    "
        write(1,fmt="(1x,a)",advance = "no") "#x"
        write(1,fmt="(1x,a)",advance = "no") "      "
        write(1,fmt="(1x,a)",advance = "no") "#y"
        write(1,fmt="(1x,a)",advance = "no") "      "
        write(1,fmt="(1x,a)",advance = "no") "#Re(psi)"
        write(1,fmt="(1x,a)",advance = "no") "  "
        write(1,fmt="(1x,a)",advance = "no") "#Im(psi)"
        write(1,fmt="(1x,a)",advance = "no") "  "
        write(1,fmt="(1x,a)",advance = "no") "#abs(psi)^2"
        write(1,*) " "

        do jj = 1,size(x_array)
                do j = 1,size(y_array)
                            write(1,FMT,advance = "no") x_array(jj)
                            write(1,fmt="(1x,a)",advance = "no") " "
                            write(1,FMT,advance = "no") y_array(j)
                            write(1,fmt="(1x,a)",advance = "no") " "
                            write(1,FMT,advance = "no") REALPART(psi_matrix(jj,j)) 
                            write(1,fmt="(1x,a)",advance = "no") " "
                            write(1,FMT,advance = "no") IMAGPART(psi_matrix(jj,j)) 
                            write(1,fmt="(1x,a)",advance = "no") " "
                            write(1,FMT,advance = "no") abs(psi_matrix(jj,j))**2 
                            write(1,*) " "
                
                end do 
            
        end do 



    end subroutine

    subroutine write_out_3index_for_gif(psi_matrix,x_array,y_array,every,delta_t)
        complex * 16,intent(in) :: psi_matrix(:,:,:)
        real * 8 ,intent(in) :: x_array(:),y_array(:) , delta_t 
        integer, intent(in):: every 
        complex * 16 :: current_psi(size(psi_matrix,2),size(psi_matrix,3))
        integer :: time_index , x_index , y_index ,shape_array(3)
        CHARACTER(LEN=20) :: FMT = "(F20.12)"
        real * 8,allocatable :: temp_info_array(:)
        integer ::max_ind

        
        max_ind = MAX(size(x_array),size(y_array))
        allocate(temp_info_array(max_ind))
        shape_array = shape(psi_matrix)
        print*,"shape array",shape_array
        open (1, file = "fourierdata_2d_gif_abs_test.dat")
        temp_info_array = 0.0_dp 

        temp_info_array(1) = REAL(size(x_array))
        temp_info_array(2) = REAL(size(y_array))
        temp_info_array(3) = REAL(every)
        temp_info_array(4) = delta_t
        temp_info_array(5) = shape_array(1)
        write (1,*) "# SIZE OF X        SIZE OF Y         EVERY               DELTA_T           NUM TIMESTEPS"
        do x_index = 1,max_ind
            write(1,FMT,advance = "no") temp_info_array(x_index)
        end do 
        write(1,*) " "
        write(1,*) "#x_array"

        do x_index = 1,max_ind
            if (x_index > size(x_array)) then 
                write(1,FMT,advance = "no") 0.0_dp
            else 
                write(1,FMT,advance = "no") x_array(x_index)
            end if
        end do 
        write(1,*) " "

        write(1,*) "#y_array"
        do y_index = 1,max_ind
            if (y_index > size(y_array)) then 
                write(1,FMT,advance = "no") 0.0_dp
            else 
                write(1,FMT,advance = "no") y_array(y_index)
            end if 
        end do 
        write(1,*) " "

        write(1,*) "####"
        do time_index = 1,shape_array(1)
            current_psi = psi_matrix(time_index,:,:)
            !print*,current_psi
            !print*,"###"
            do y_index = 1,size(y_array)
                do x_index = 1,size(x_array)
                    write(1,FMT,advance = "no") abs(current_psi(x_index,y_index))
                end do 
                write(1,*) " "
                
            end do 
            write(1,*) "#######################"
        end do 


    end subroutine

    function gaussian_2d(x_array,y_array,width,wavenumber) result(gaussian_2d_array)
        real * 8 , intent(in) :: x_array(:),y_array(:),width,wavenumber
        complex * 16 :: gaussian_2d_array(size(x_array),size(y_array))
        integer :: i,j 

        do i = 1,size(x_array)
            do j = 1,size(y_array)
                !fortran is column majored remember
                gaussian_2d_array(i,j) = EXP(-(x_array(i)**2 + y_array(j)**2)/width)*EXP(complex(0.0_dp,x_array(i)*wavenumber))
            end do 
        end do


    end function 

    function gaussian_2d_potential(x_array,y_array,width,centerx,centery) result(gaussian_2d_array)
        real * 8 , intent(in) :: x_array(:),y_array(:),width,centerx,centery
        real * 8 :: gaussian_2d_array(size(x_array),size(y_array))
        integer :: i,j 

        do i = 1,size(x_array)
            do j = 1,size(y_array)
                !fortran is column majored remember
                gaussian_2d_array(i,j) = EXP(-((x_array(i)-centerx)**2 + (y_array(j)-centery)**2)/width)
            end do 
        end do


    end function 


end module tdse_fourier_2d