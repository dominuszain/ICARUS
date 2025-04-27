module icarus_utils
    
    implicit none
    
    real :: pi = 3.14159265358979
    real :: e = 2.71828182845904
    real :: c = 2.997e8
    real :: avo = 6.022e23
    real :: mass_neutron = 939.6
    real :: amu = 931.5
    
    contains
    
        subroutine verify_input(col1, col2, col3, flag)
        
            implicit none
            
            logical, intent(out) :: flag
            real, dimension(:), intent(in) :: col2
            character(len = 12), dimension(:), intent(in) :: col1, col3
            
            if (.not.(col1(1) == "target" .and. col1(2) == "points" .and. col1(3) == "step" .and. col1(4) == "path")) then
            
                print *, "input file error: incorrect format."
                print *, "problem lies somewhere in column 1."
                flag = .false.
            
            else if (.not.(col3(1) == 'n' .and. col3(2) == 'n' .and. col3(3) == 'n')) then
            
                print *, "input file error: incorrect format."
                print *, "problem lies somewhere in column 3"
                flag = .false.
            
            else if (.not.(col2(4) == -1)) then
            
                print *, "input file error: incorrect format."
                print *, "problem lies somewhere in column 2."
                flag = .false.
            
            else
            
                flag = .true.
                
            end if
        
        end subroutine

end module icarus_utils


program icarus
    
    use icarus_utils
    
    implicit none
    
    logical :: flag = .false.
    real :: target, red_m, step, v_t, const, mxw, rr
    integer :: i, j, start_e, end_e, step_e, col_n, max_ener
    character(len = 12), dimension(:), allocatable :: col1, col3
    real, dimension(:), allocatable :: col2, output
    real, dimension(:, :), allocatable :: input
    character(len = 12) :: inp_name
    
    col_n = 4
    
    start_e = 5
    end_e = 120
    step_e = 5
    
    allocate(col1(col_n), col2(col_n), col3(col_n))
    
    open(1, file = "icarus.inp")
    
    do i = 1, col_n
    
        read (1, *) col1(i), col2(i), col3(i)
    
    end do
    
    close(1)
    
    call verify_input(col1, col2, col3, flag)
    
    if (flag .eqv. .false.) then
    
        stop
    
    end if
    
    inp_name = col3(4)
    max_ener = int(col2(2))
    step = col2(3)
    
    target = col2(1)
    
    red_m = (mass_neutron * target * amu) / (((target) * (amu)) + mass_neutron)
    
    allocate(input(max_ener, 2), output(max_ener))
    
    open(1, file = inp_name)
    
    do i = 1, max_ener
    
        read (1, *) input(i, :)
    
    end do
    
    close(1)
    
    open(1, file = "icarus.out")
    
    write (1, *) "ICARUS is an s- process nucleosynthesis code that calculates MACS and RR."
    write (1, *) "Authored by Zain Ul Abideen (dominuszain@gmail.com) as part of the FYP."
    write (1, *) ""
    write (1, *) "col1: Energy (keV)"
    write (1, *) "col2: MACS (mb)"
    write (1, *) "col3: RR (cm3 s-1 mol-1)"
    write (1, *) ""
    
    do i = start_e, end_e, step_e

        v_t = sqrt(2 * (real(i) / 1000) * c**2 / red_m)
        const = 2 / (sqrt(pi) * (real(i) / 1000)**2)
        
        do j = 1, max_ener
        
            output(j) = const * input(j, 2) * input(j, 1) * e**(-input(j, 1) / (real(i) / 1000)) * step
        
        end do
        
        mxw = sum(output)
        rr = mxw * v_t * avo * 1e-25
        
        write (1, *) i, mxw, rr
    
    end do
    
    write (1, *)
    print *, "ICARUS congratulates you on performing a successful calculation."
    write (1, *) "ICARUS congratulates you on performing a successful calculation."
    
    close(1)
    
    deallocate(input, output, col1, col2, col3)

end program icarus
