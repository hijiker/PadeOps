    subroutine test_lattice_definition(sout)
        real(rkind) :: s
        integer :: idx
        real(rkind), intent(out) :: sout

        sout = zero

        call message(0,"PERFORMING LATTICE TESTS")
        ! Test 1: sum of weights
        s = sum(w); call message(1,"TEST 1: Sum of weights ",s)
        sout = sout + s - one

        ! Test 2: wi ci
        s = sum(w*cx); call message(1,"TEST 2a: Sum of {w_i c_i_x} ",s); sout = sout + s
        s = sum(w*cy); call message(1,"TEST 2b: Sum of {w_i c_i_y} ",s); sout = sout + s
        s = sum(w*cz); call message(1,"TEST 2c: Sum of {w_i c_i_z} ",s); sout = sout + s

        ! Test 3: wi ci_alpha ci_beta
        s = sum(w*cx*cx); call message(1,"TEST 3a: Sum of {w_i c_i_x c_i_x} ",s); sout = sout + s - csq
        s = sum(w*cx*cy); call message(1,"TEST 3a: Sum of {w_i c_i_x c_i_y} ",s); sout = sout + s
        s = sum(w*cx*cz); call message(1,"TEST 3a: Sum of {w_i c_i_x c_i_z} ",s); sout = sout + s
        s = sum(w*cy*cy); call message(1,"TEST 3a: Sum of {w_i c_i_y c_i_y} ",s); sout = sout + s - csq
        s = sum(w*cy*cz); call message(1,"TEST 3a: Sum of {w_i c_i_y c_i_z} ",s); sout = sout + s
        s = sum(w*cz*cz); call message(1,"TEST 3a: Sum of {w_i c_i_z c_i_z} ",s); sout = sout + s - csq

        ! Test 4: Boundary identification
        s= zero 
        do idx = 1,size(fparallel)
            s = s + cz(fparallel(idx))
        end do 
        call message(1,"TEST 4: Sum fparallel (should be zero):", s)
        sout = sout + s

        s= zero 
        do idx = 1,size(fplus)
            s = s + cz(fplus(idx)) - one
        end do 
        call message(1,"TEST 4: Sum fplus (should be zero):", s)
        sout = sout + s


        s= zero 
        do idx = 1,size(fminus)
            s = s + cz(fminus(idx)) + one
        end do 
        call message(1,"TEST 4: Sum fminus (should be zero):", s)
        sout = sout + s

    
        s= zero 
        do idx = 1,size(fminus)
            s = s + cz(fminus(idx)) + cz(fplus(idx))
        end do 
        call message(1,"TEST 4: Sum foppo (should be zero):", s)
        sout = sout + s



    end subroutine 
