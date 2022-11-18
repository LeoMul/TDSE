program main 
    real * 8 :: psi_array(3,2,2)
    real * 8 ::edit_array(2,2)
    psi_array = 1.0
    edit_array = 2.0

    psi_array(1,:,:) = edit_array

    print*,psi_array(1,:,:)

end program 