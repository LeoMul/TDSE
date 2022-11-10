program main
    use tdse

    implicit none 
    integer,parameter:: n = 3
    complex*16 :: m_matrix(n,n),v_array(n),alpha,delta_t ,sol(n),test(n)
    integer :: i

    !not necessary
    delta_t = (1.0_dp,0.0_dp)
    alpha = (1.0_dp,0.0_dp)
    v_array = (1.0_dp,0.0_dp)
    m_matrix = create_m_matrix(v_array,alpha,delta_t)


    !check
    m_matrix = reshape((/(1.0_dp,0.0_dp),(1.0_dp,0.0_dp),(0.0_dp,0.0_dp),(20.0_dp,0.0_dp),(1.0_dp,0.0_dp),(2.0_dp,0.0_dp),(0.0_dp,0.0_dp),(1.0_dp,0.0_dp),(1.0_dp,0.0_dp)/),(/n,n/))

    print*,"begin matrix"
    do i = 1,n
        print*, REAL(m_matrix(i,:))
    end do
    print*,"end of matrix"

    print*,"rhs"
    v_array = (/(1.0_dp,0.0_dp),(3.0_dp,0.0_dp),(9.0_dp,0.0_dp)/)
    print*,REAL(v_array)

    sol = fast_tri_diag_solution(m_matrix,v_array)
    test = fast_tria_multiply_by_vector(m_matrix,sol)
    print*,"sol:"
    print*, REAL(sol)
    print*,"check:"
    print*,REAL(test)


end program main