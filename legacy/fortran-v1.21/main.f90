module icarus_utils
    
    implicit none
    
    real, parameter :: pi = 3.14159265358979
    real, parameter :: e = 2.71828182845904
    real, parameter :: c = 2.997e8
    real, parameter :: avo = 6.022e23
    real, parameter :: mass_neutron = 939.6
    real, parameter :: amu = 931.5
    
    contains
    
        subroutine cool_logo()
        
        implicit none
        
        print *,"_________ _______  _______  _______           _______ "
        print *,"\__   __/(  ____ \(  ___  )(  ____ )|\     /|(  ____ \"
        print *,"   ) (   | (    \/| (   ) || (    )|| )   ( || (    \/"
        print *,"   | |   | |      | (___) || (____)|| |   | || (_____ "
        print *,"   | |   | |      |  ___  ||     __)| |   | |(_____  )"
        print *,"   | |   | |      | (   ) || (\ (   | |   | |      ) |"
        print *,"___) (___| (____/\| )   ( || ) \ \__| (___) |/\____) |"
        print *,"\_______/(_______/|/     \||/   \__/(_______)\_______) v1.2"
        print *, ""
        print *, "Authored by Zain Ul Abideen as part of FYP (BS-SS09-IST)"
        print *, "contact: dominuszain@gmail.com"
        
        
        end subroutine

end module icarus_utils


program icarus
    
    use icarus_utils
    
    implicit none
    
    real :: tgt, red_m, step, v_t, const, mxw, rr
    integer :: i, j, start_e, end_e, step_e, max_ener
    real, dimension(:), allocatable :: output
    real, dimension(:, :), allocatable :: input
    character(len=32) :: inp_name
    
    start_e = 5
    end_e = 100
    step_e = 5
    
    call cool_logo()
    
    read *, tgt, max_ener
    read *, inp_name
    
    red_m = (mass_neutron * tgt * amu) / (((tgt) * (amu)) + mass_neutron)
    
    allocate(input(max_ener, 2), output(max_ener))
    
    open(1, file = inp_name)
    
    do i = 1, max_ener
    
        read (1, *) input(i, :)
    
    end do
    
    close(1)
    
    print *, ""
    
    do i = start_e, end_e, step_e

        v_t = sqrt(2 * (real(i) / 1000) * c**2 / red_m)
        const = 2 / (sqrt(pi) * (real(i) / 1000)**2)
        
        do j = 1, max_ener - 1
        	
        	step = input(j + 1, 1) - input(j, 1)
            output(j) = const * input(j, 2) * input(j, 1) * e**(-input(j, 1) / (real(i) / 1000)) * step
        
        end do
        
        mxw = sum(output)
        rr = mxw * v_t * avo * 1e-25
        
        print *, i, 1.1605*i/100, mxw, rr
    
    end do
    
	print *, ""
    print *, "ICARUS congratulates you on performing a successful calculation."
    
    
    deallocate(input, output)

end program icarus
