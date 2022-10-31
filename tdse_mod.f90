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

    function gaussian(x_array,wavenumber,delta_x) result (gauss)
        real*8,intent(in)::x_array(:),wavenumber,delta_x
        complex*16 :: gauss(size(x_array))
        integer:: j 
        complex*16:: g
        do j = 1,size(gauss)
            g = EXP(-x_array(j)*x_array(j)/0.2_dp+im_unit*wavenumber*x_array(j))
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

end module tdse 