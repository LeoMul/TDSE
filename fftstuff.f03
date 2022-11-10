module FFT_MOD
    use, intrinsic :: iso_c_binding
    include 'fftw3.f03'
    contains 

    function fft(in) result(out)
        complex(c_double_complex),intent(inout) :: in(:)
        complex(c_double_complex) :: out(size(in))
        type(C_PTR) :: p
        p = fftw_plan_dft_1d(size(in), in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        call fftw_execute_dft(p, in, out)
        call fftw_destroy_plan(p)


    end function

    function ifft(in) result(out)

        complex(c_double_complex),intent(inout) :: in(:)
        complex(c_double_complex) :: out(size(in))
        type(C_PTR) :: p
        p = fftw_plan_dft_1d(size(in), in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
        call fftw_execute_dft(p, in, out)
        out = out / size(out)
        call fftw_destroy_plan(p)

    end function
    
    function create_fftfreq(n,d) result (fftfreq)
        !makes the correct k's (needs multiplied by 2pi) and scales by length.

        real * 8,intent(in) :: d 
        integer, intent(in) :: n 
        real * 8 :: fftfreq(n)
        integer :: i 

        if (modulo(n,2) == 0) then
            do i = 1,n/2 
                fftfreq(i) = i-1
            end do 
            do i = n,n/2+1,-1 
                fftfreq(i) = i-1-n
            end do

        else 
            do i = 1,(n+1)/2 
                fftfreq(i) = i-1
            end do 
            do i = n,(n+1)/2+1,-1 
                fftfreq(i) = i-1-n
            end do
        end if 
        fftfreq = fftfreq / (n*d)
    end function create_fftfreq



end module