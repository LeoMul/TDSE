program fourier_test
    use tdse
    use FFT_MOD
    use, intrinsic :: iso_c_binding
    implicit none
    CHARACTER(LEN=20) :: FMT = "(F20.12)"
    real*8 :: pi = 3.14159265359
    integer, parameter :: N = 1000
    integer :: j,l
    real*8::x_array(N),x_0, x_N,k_array(N),delta_x
    complex * 16 ::  y_array(N)
    complex * 16 :: test_array(N),sum
    complex(c_double_complex),dimension(N) :: out,deriv

    x_0 = -10.0_dp 
    x_N = 10.0_dp
    x_array = my_linspace(x_0,x_N,N)
    delta_x = (x_array(2)-x_array(1))

    k_array = 2.0_dp*pi*create_fftfreq(N,delta_x)
    y_array = gaussian(x_array,0.0_dp,0.1_dp)

    

    open (1, file = "fftwtest.dat")
    
    


    out = -k_array*k_array*fft(y_array)
    deriv = ifft(out)

    y_array = y_array*(4.0_dp*x_array*x_array-2.0_dp)

    do j = 1,N 
        sum = COMPLEX(0.0_dp,0.0_dp)
        do l = 1,N
            sum = sum + exp(im_unit*(j-1)*(l-1)*2.0*pi/N) * out(l)
        end do 
        test_array(j) = sum
    end do 
    test_array = test_array/N
    
    do j = 1,size(x_array)

        write(1,FMT,advance = "no") x_array(j)
        write(1,fmt="(1x,a)",advance = "no") " "
        write(1,FMT,advance = "no") REAL(y_array(j))
        write(1,fmt="(1x,a)",advance = "no") " "
        write(1,FMT,advance = "no") REALPART(deriv(j))
        write(1,*) " "
         
    end do 


end program fourier_test