module tdse
    
    integer, parameter :: dp = selected_real_kind(15, 307)
    complex*16 :: im_unit = (0.0_dp,1.0_dp)

    contains 

    function calculate_psi_next_time_step(psi_array_in,left_boundary,right_boundary,V_array,alpha,delta_t,delta_x) result(psi_array_out)
        complex*16,intent(in) :: psi_array_in(:),left_boundary,right_boundary
        integer :: n ,ii
        complex*16 :: psi_array_out(size(psi_array_in))
        real*8,intent(in)::alpha,V_array(:),delta_t,delta_x
        real*8 :: norm

        n = size(psi_array_in)
        psi_array_out(1) = left_boundary
        psi_array_out(n) = right_boundary
        do ii = 2,n-1
            psi_array_out(ii) =  psi_array_in(ii) - im_unit * delta_t * V_array(ii)*psi_array_in(ii) + 0.5_dp * im_unit * alpha * (psi_array_in(ii+1)+psi_array_in(ii-1)-2.0_dp*psi_array_in(ii))
            !print*, 0.5_dp * im_unit * alpha * (psi_array_in(ii+1)+psi_array_in(ii-1)-2.0_dp*psi_array_in(ii))
        end do 
        
        norm  = sqrt(sum(abs(psi_array_out)**2) * delta_x)

        psi_array_out = psi_array_out/norm
        
    end function

    function gaussian(x_array,wavenumber,width,offset) result (gauss)
        real*8,intent(in)::x_array(:),wavenumber,width,offset
        real *8 ::new_x_array(size(x_array))
        complex*16 :: gauss(size(x_array))
        integer:: j 
        complex*16:: g
        new_x_array = x_array - offset
        do j = 1,size(gauss)
            g = COMPLEX(EXP((-new_x_array(j)*new_x_array(j)/width)),0.0_dp)*EXP(im_unit*wavenumber*x_array(j))
            gauss(j) = g
        end do
        !gauss = gauss/ sqrt((sum(abs(gauss)**2) * delta_x))
    end function gaussian
    

    function create_v_array(x_array) result(v)
        !free particle for now
        real*8, intent(in) :: x_array(:)
        real*8 :: v(size(x_array))
        v = 0.0_dp

    end function 


    function my_linspace(x_0,x_last,N)
        !makes an evenly spaced array of real numbers betweenn x_0 and x_last of length N
        real*8, intent(in)::x_0,x_last
        integer::N
        real*8::h
        real*8::my_linspace(N)

        h = (x_last-x_0)/(N-1)

        my_linspace = my_arange(x_0,x_last,h)

    end function my_linspace

    function my_arange(x_0,x_last,h) result (rangearray)
        
        real*8, intent(in)::x_0,x_last,h
        integer::N,i
        real*8,allocatable::rangearray(:)
        
        N = NINT((x_last-x_0)/h)+ 1
        allocate(rangearray(N))
        
        do i = 1,N
            rangearray(i) = x_0 + (i-1)*h
        end do

    end function my_arange

    function create_m_matrix(v_array_without_ends,alpha,delta_t) result(m_matrix)
        complex*16,intent(in) :: v_array_without_ends(:)
        complex*16 :: m_matrix(size(v_array_without_ends),size(v_array_without_ends))
        complex*16,intent(in) :: alpha ,delta_t

        integer :: k,j,n 

        n = size(v_array_without_ends)

        do k = 1,n 
            do j = 1,n 
                if (j == k) then 
                    m_matrix(j,k) = im_unit * (alpha +0.5_dp*delta_t*v_array_without_ends(j))
                
                else if ((j == k-1) .or. (j == k+1))  then 
                    m_matrix(j,k) = -0.5_dp*im_unit*alpha
                
                else
                    m_matrix(j,k) = (0.0_dp,0.0_dp)
                end if 

            end do 
        end do
    end function create_m_matrix

    function fast_tri_diag_solution(tri_matrix,d) result(solution)
        !https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm

        !O(n) solution to tria_diag * x = b
        complex*16, intent(in) :: d(:) !rhs
        complex*16,intent(in):: tri_matrix(:,:)
        complex*16 :: solution(size(d)),a,b,c

        complex*16 :: d_prime (size(d)),c_prime(size(d)-1)

        integer :: j,n 

        n = size(d)
        
        b = tri_matrix(1,1)

        d_prime(1) = d(1)/b
        c_prime(1) = tri_matrix(1,2)/b

        do j = 2,n-1
            a = tri_matrix(j,j-1)
            b = tri_matrix(j,j)

            c = tri_matrix(j,j+1)

            c_prime(j) = c/ (b-a*c_prime(j-1))

            d_prime(j) = (d(j) - a*d_prime(j-1))/(b-a*c_prime(j-1))
        end do 
        !to save having an if in the last one
        a = tri_matrix(n,n-1)
        b = tri_matrix(n,n)
        d_prime(n) = (d(n) - a*d_prime(n-1))/(b-a*c_prime(n-1))


        solution(n) = d_prime(n)
        do j = n-1,1,-1
            solution(j) = d_prime(j) - c_prime(j)* solution(j+1)
        end do 


    end function fast_tri_diag_solution

    function fast_tria_multiply_by_vector(tri_matrix,a) result (b)
        complex * 16 ,intent(in) :: tri_matrix(:,:),a(:)
        complex * 16 :: b(size(a))
        integer :: j , n 

        n = size(a)

        b(1) = tri_matrix(1,1)*a(1) + tri_matrix(1,2)*a(2)

        b(n) = tri_matrix(n,n-1)*a(n-1)+tri_matrix(n,n)*a(n) 

        do j = 2,n-1
            b(j) = tri_matrix(j,j-1)*a(j-1)+ tri_matrix(j,j)*a(j) + tri_matrix(j,j+1)*a(j+1)
        end do 

        

    end function

    function create_indentity_matrix(n) result (ident)
        integer,intent(in) :: n 
        complex*16 :: ident(n,n) 
        integer :: j,k 

        do j = 1,n
            do k = 1,n 
                if  (j == k) then 
                    ident(j,j) = (1.0_dp,0.0_dp)
                else 
                    ident(j,k) = (0.0_dp,0.0_dp)

                end if 
            end do  
        end do 
        

    end function create_indentity_matrix

    function schrodinger_crank_time_independent_one_step(psi_array,m_matrix,left_bound,right_bound,identity,alpha,delta_x) result (new_psi)
        complex * 16 ,intent(in):: psi_array(:) ,left_bound,right_bound,identity(:,:),m_matrix(:,:),alpha
        complex * 16 :: left_matrix(size(psi_array)-2,size(psi_array)-2), right_matrix(size(psi_array)-2,size(psi_array)-2),rhs(size(psi_array)-2)
        complex * 16 :: new_psi (size(psi_array)),boundary(size(psi_array)-2)
        real * 8 ,intent(in) :: delta_x
        real * 8 :: norm_factor
        integer :: n 

        n = size(psi_array)
        !COMPLEX(alpha,0.0_dp),COMPLEX(delta_t,0.0_dp)
        new_psi(1) = left_bound
        new_psi(n) = right_bound
        boundary = 0.0_dp 
        boundary(1) = left_bound
        boundary(n-2) = right_bound
        !check this correctly incorporates the boundaries
        boundary = alpha*boundary
        !technically the boundary condit

        right_matrix = identity - m_matrix
        left_matrix = identity + m_matrix
        rhs = fast_tria_multiply_by_vector(right_matrix,psi_array(2:n-1)) + boundary


        new_psi(2:n-1) = fast_tri_diag_solution(left_matrix,rhs)
        norm_factor = sum(abs(new_psi)**2)*delta_x
        new_psi = new_psi / sqrt(norm_factor)

    end function schrodinger_crank_time_independent_one_step

    function create_square_barrier(x_array,start,length,strength) result(v_array)
        real * 8,intent(in):: x_array(:),start,length,strength
        complex*16 :: v_array(size(x_array)) 
        integer :: j 

        do j = 1,size(x_array)
            if ((x_array(j) .ge. start) .and. (x_array(j) .le. start+length) )then
                v_array(j) = COMPLEX(strength , 0.0_dp)
            else  
                v_array(j) = COMPLEX(0.0_dp,0.0_dp)
            end if
        end do 

    end function create_square_barrier

    function nsd(delta_x,psi_array) result (nsd_array)
        real * 8,intent(in) :: delta_x 
        complex * 16 , intent(in) :: psi_array (:)
        complex * 16 :: nsd_array (size(psi_array)-2)
        integer :: j 

        do j = 1,size(psi_array)-2
            nsd_array(j) = psi_array(j) + psi_array(j+2) - COMPLEX(2.0_dp,0.0_dp) * psi_array(j+1)
        end do 
        nsd_array = nsd_array/COMPLEX(delta_x**2,0.0_dp)
    end function nsd

    function calculate_energy_integral(delta_x,v_array,psi_array) result(integral)
        complex * 16, intent(in) :: v_array(:),psi_array(:)
        real * 8, intent(in) :: delta_x
        complex * 16 :: integral 
        integer :: n 

        n = size(psi_array)
        !print*,nsd(delta_x,psi_array)
        integral = - complex(0.5_dp,0.0_dp)* DOT_PRODUCT(nsd(delta_x,psi_array),conjg(psi_array(2:n-1)))
        integral = integral + dot_product(v_array(2:n-1),abs(psi_array(2:n-1))**2)
        !print*,integral
        integral = delta_x * integral / (sqrt(sum(abs(psi_array(2:n-1))**2))*delta_x)




    end function calculate_energy_integral

    

    

    subroutine rearrange_k(k_array) 
        real * 8, intent(inout) :: k_array(:)
        integer :: n 
        real * 8 :: temp_ks(size(k_array))

        n = size(k_array)

        temp_ks = k_array 

        k_array(1:n/2 ) = temp_ks ((n/2 + 1): n)
        k_array((n/2 + 1 ):n) = temp_ks(1:n/2 )

        temp_ks = k_array

        !k_array(1) = temp_ks(n)
        !k_array(2:n) = temp_ks(1:n-1)

    end subroutine rearrange_k

    subroutine rearrange_fft(k_array) 
        complex * 8, intent(inout) :: k_array(:)
        integer :: n 
        complex * 8 :: temp_ks(size(k_array))

        n = size(k_array)

        temp_ks = k_array 

        k_array(1:n/2 ) = temp_ks ((n/2 + 1): n)
        k_array((n/2 + 1 ):n) = temp_ks(1:n/2 )

        !temp_ks = k_array
!
        !k_array(1) = temp_ks(n)
        !k_array(2:n) = temp_ks(1:n-1)

    end subroutine rearrange_fft







end module tdse 