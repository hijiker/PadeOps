module NS_1d_rhs
    use kind_parameters, only: rkind
    use cd06staggstuff, only:cd06stagg
    use interpolation
    implicit none

    real(rkind), parameter :: dz_match = 0.10d0 !50.d0/1000.d0
    integer, parameter :: nz = 40
    real(rkind), dimension(:), allocatable :: nu_t, z
    real(rkind) :: umatch = 0.d0 !14.8439612264d0

    real(rkind) :: ustar = 1.d0 
    integer :: scheme = 2 ! 1: 2nd order, 2: CD06
    real(rkind), parameter :: Retau = 5000
    
    type(cd06stagg) :: der
    real(rkind) :: dz
    real(rkind), dimension(:), allocatable :: bsp, csp, dsp 
    real(rkind), dimension(:,:), allocatable :: data2read
    
    real(rkind), dimension(:), allocatable :: zplus, yplus

contains
    
    subroutine get_nu_t(ustar)
        real(rkind), intent(in) :: ustar
        integer :: idx, nzdat

        zplus = (z + 1)*Retau*ustar
        
        nzdat = size(data2read,1)
        !call spline (yplus, data2read(:,2), bsp, csp, dsp, size(bsp))
        do idx = 1,nz/2
            if (zplus(idx)>Retau) then
                nu_t(idx) = data2read(nzdat/2,2)
            else
                nu_t(idx) = ispline(zplus(idx), yplus, data2read(:,2), bsp, csp, dsp, size(bsp)) 
            end if 
        end do 
            
        do idx = nz/2+1,nz
            nu_t(idx) = nu_t(nz-idx+1)
        end do 

    end subroutine 


    subroutine get_RHS(nz, u, rhs, dudz, F, nu_t, tmpC1, tmpC2, tmpE)
        use constants, only: zero 
        integer, intent(in) :: nz
        real(rkind), intent(in) :: F
        real(rkind), intent(in) , dimension(1,1,nz+1) :: u
        real(rkind), intent(out), dimension(1,1,nz+1) :: rhs, dudz, tmpE
        real(rkind), intent(out), dimension(1,1,nz) :: tmpC1, tmpC2
        real(rkind), intent(in), dimension(nz) :: nu_t
        integer :: idx
        real(rkind) :: ustarnew, nu_wall

        rhs = 0.d0

        ! STEP 1: Get ustar
        ustarnew = find_utau(ustar, u(1,1,2), dz)
        ustar = ustarnew 

        ! STEP 2: Get dudz
        call der%InterpZ_E2C(u,tmpC1,1,1)
        call der%ddz_C2E(tmpC1,dudz,1,1)
        
        ! STEP 3: Get nu_wall
        nu_wall = ustar*ustar/(dudz(1,1,1)+1.D-14)
       
        ! STEP 4: Compute total viscosity and replace boundary
        do idx = 1,nz+1
            tmpE(1,1,idx) = nu_t(idx) + (1.d0/Retau)
        end do 
        tmpE(1,1,1) = nu_wall
        tmpE(1,1,nz+1) = nu_wall

        ! STEP 5: Interpolate viscosity
        call der%InterpZ_E2C(tmpE,tmpC1,1,1)
        
        ! STEP 6: Interpolate derivative
        call der%ddz_E2C(u,tmpC2,1,1)

        ! STEP 7: Compute stress
        tmpC2 = tmpC2*tmpC1

        ! STEP 8: Compute RHS
        call der%ddz_C2E(tmpC2,rhs,1,1)

        ! STEP 9: Add body force
        rhs = rhs + F

        ! STEP 10: 
        rhs(1,1,1) = zero
        rhs(1,1,nz+1) = zero

    end subroutine 

    pure elemental function get_musker_profile(yp) result(up)
        real(rkind), intent(in) :: yp
        real(rkind) :: up

        up = 5.424d0*atan(0.11976047904*yp - 0.48802395209580833) + 0.434d0*log(((yp + 10.6)**9.6)/((yp**2 - 8.15*yp + 86)**2)) - 3.50727901936264842

    end function

    function find_utau(utau_guess, umatch, zmatch) result(utau)
        real(rkind), intent(in) :: utau_guess, umatch, zmatch
        real(rkind) :: zplus, utau_diff, up1, utaunew
        real(rkind) :: utau

        utau_diff = utau_guess
        utau = utau_guess
        do while (utau_diff > 1.d-8)
            zplus = zmatch*Retau*utau 
            up1 = get_musker_profile(zplus)                
            utaunew = umatch/up1
            utau_diff = abs(utaunew - utau)
            utau = utaunew
        end do 

        ! print*, tau
    end function 

        
    function get_utarget(u1,z1,z0) result(u0)
        real(rkind), intent(in) :: u1, z1, z0
        real(rkind) ::  up0
        real(rkind) :: u0, ustarnew

        ustarnew = find_utau(ustar, u1, z1)
        ustar = ustarnew

        up0 = get_musker_profile(z0*Retau*ustar)
        u0 = up0*ustar

    end function 
end module 


program FV_testing_RANS
    use kind_parameters, only: rkind, clen 
    use basic_io
    use NS_1d_rhs
    use gridtools, only: linspace

    implicit none

    real(rkind), parameter :: tstop = 5000.0d0 
    real(rkind), parameter :: cvisc =  0.4d0
    real(rkind), dimension(:,:,:), allocatable :: dt 
    real(rkind), dimension(:,:,:), allocatable :: u, tmpC1, rhs, dudz,  tmpC2, tmpE
    real(rkind), dimension(:,:,:), allocatable :: uin, k1, k2, k3, k4
    real(rkind), parameter :: F = 1.d0 
    character(len=clen) :: fname, input_fname = "/scratch/04076/tg833754/Cebecci_smith/Cebeci_smith_ReTau5000.txt"
    integer :: idx, nzdat, tprint = 100000, tdatadump = 5000000 
    integer, parameter :: tscheme = 3
    real(rkind), dimension(:,:), allocatable :: data2read2

    ! get the nu_t from disk
    call read_2d_ascii(data2read,input_fname)

    allocate(bsp(size(data2read,1)))
    allocate(csp(size(data2read,1)))
    allocate(dsp(size(data2read,1)))
    allocate(yplus(size(data2read,1)))


    allocate(z(nz+1),nu_t(nz+1))
    allocate(u(1,1,nz+1))
    allocate(dudz(1,1,nz+1))
    allocate(tmpC1(1,1,nz))
    allocate(tmpC2(1,1,nz))
    allocate(rhs(1,1,nz+1))
    allocate(zplus(nz+1))
    allocate(tmpE(1,1,nz+1))

    allocate(dt(1,1,nz+1))
   
    z = linspace(-1.d0,1.d0,nz+1)
    dz = z(2) - z(1)
    call get_nu_t(ustar)

    call der%init(nz, dz, isTopEven=.false., isBotEven=.false., isTopSided=.true., isBotSided=.true.)
    
    u(1,1,:) = 1.d0*(1 - z**2)  !am 

    yplus = (data2read(:,1)+1)*Retau
    zplus = (z + 1)*Retau*ustar
    call spline (yplus, data2read(:,2), bsp, csp, dsp, size(bsp))
    nzdat = size(data2read,1)
    do idx = 1,nz+1
           nu_t(idx) = ispline(zplus(idx), yplus, data2read(:,2), bsp, csp, dsp, size(bsp)) 
    end do 
    print*, "Nz:", nz

    do idx = 1,nz
        dt(1,1,idx) = cvisc*dz*dz/(1.d0/Retau + nu_t(idx))
    end do 
    dt = 1.d-4

    allocate(data2read2(nz+1,2))

    allocate(k1(1,1,nz+1))
    allocate(k2(1,1,nz+1))
    allocate(k3(1,1,nz+1))
    allocate(k4(1,1,nz+1))
    allocate(uin(1,1,nz+1))

    deallocate(data2read)
    allocate(data2read(nz+1,2))
    idx = 1
    do while(idx<5000000)
        ! Interpolate cell to edge

        select case (tscheme)
        case (0)
            call get_RHS(nz, u, k1, dudz, F, nu_t, tmpC1, tmpC2, tmpE)
            u = u + k1*dt
        case (1)
            uin = u
            call get_RHS(nz, uin, k1, dudz, F, nu_t, tmpC1, tmpC2, tmpE)
            k1 = dt*k1

            uin = u + 0.5d0*k1
            call get_RHS(nz, uin, k2, dudz, F, nu_t, tmpC1, tmpC2, tmpE)
            k2 = dt*k2

            uin = u + 0.5d0*k2
            call get_RHS(nz, uin, k3, dudz, F, nu_t, tmpC1, tmpC2, tmpE)
            k3 = dt*k3

            uin = u + k3
            call get_RHS(nz, uin, k4, dudz, F, nu_t, tmpC1, tmpC2, tmpE)
            k4 = dt*k4

            u = u + (1.d0/6.d0)*(k1 + k2 + k3 + k4)
        case (2)
            call get_RHS(nz, u, k1, dudz, F, nu_t, tmpC1, tmpC2, tmpE)
            k1 = u + dt*k1

            call get_RHS(nz, k1, k2, dudz, F, nu_t, tmpC1, tmpC2, tmpE)
            k2 = (3.d0/4.d0)*u + (1.d0/4.d0)*k1 + (1.d0/4.d0)*dt*k2

            call get_RHS(nz, k2, k3, dudz, F, nu_t, tmpC1, tmpC2, tmpE)
            
            u = (1.d0/3.d0)*u + (2.d0/3.d0)*k2 + (2.d0/3.d0)*dt*k3
        case (3)
            call get_RHS(nz, u, k1, dudz, F, nu_t, tmpC1, tmpC2, tmpE)
            k1 = u + dt*k1

            call get_RHS(nz, k1, k2, dudz, F, nu_t, tmpC1, tmpC2, tmpE)
            
            u = 0.5d0*u + 0.5d0*k1 + 0.5d0*dt*k2
        end select

        !t = t + dt
        idx = idx + 1
        if (mod(idx,tprint) == 0) then
            print*, "ubottom:", u(1,1,1)
            print*, "uquarte:", u(1,1,nz/4)
            print*, "ucenter:", u(1,1,nz/2)
            print*, "-----------" 
            print*, "utau:", ustar
            print*, "maxvalu:", maxval(u)
            print*, "minvalu:", minval(u)
            print*, "-----------" 
            print*, "-----------" 
            !print*,"Laminar error:", maxval(abs(u - ulam))
        end if

        if (mod(idx,tdatadump) == 0) then
            write(fname,"(A7,I10.10,A4)") "u_CS_ED",idx, ".dat"
            data2read2(:,1) = (z+1)*ustar*Retau
            data2read2(:,2) = u(1,1,:)/ustar
            call write_2d_ascii(data2read2,trim(fname))
        end if 
    end do 

end program 
