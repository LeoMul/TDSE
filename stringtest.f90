program stringtest
    character (len=50) :: input_trim
    real * 8 :: float_to_add 

    float_to_add = 1.0d0
    !input_trim = "qqwewqe"
    print*,float_to_add
    !input_trim = input_trim // float_to_add

    !read(input_trim,*) float_to_add
    write(input_trim,'(a7,F20.3)') "qqwewqe",float_to_add

    print*,input_trim
end program 