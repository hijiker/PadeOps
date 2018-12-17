module d3q19_testing
    use kind_parameters, only: rkind
    use constants, only: one
    implicit none

    integer :: testing_type = 0
    real(rkind) :: Reynolds_number = one/(0.03472237895076d0), Force = 2.4112862d0
    real(rkind), dimension(:,:,:), allocatable :: utrue, vtrue, wtrue
    real(rkind) :: omega = 2._rkind
    real(rkind) :: p0 = 48.225744d0
    real(rkind), dimension(:), allocatable :: z
    real(rkind), dimension(19) :: mvec_mrt_uniform 
contains
    subroutine get_womersley_exact(y,omega,nu,t, u)
        use constants, only: imi, two
        real(rkind), intent(in)  :: y, omega, nu, t
        real(rkind), intent(out) :: u
        
        real(rkind) :: alpha, D

        D = one
        alpha = (D/2)*sqrt(omega/nu)

        u = real((p0/(imi*omega))*(one-cosh(sqrt(two)*(one + imi)*alpha*y/D)/cosh(sqrt(two)*(one+imi)*alpha))*exp(imi*omega*t),rkind)

    end subroutine

    subroutine update_utrue_wormersley(time)
        real(rkind), intent(in) :: time
        integer :: idx, nz
        real(rkind) :: utmp 

        nz = size(utrue,3)
        do idx = 1,nz
            call get_womersley_exact(z(idx),omega,one/Reynolds_number,time,utmp)
            utrue(:,:,idx) =  utmp
        end do 



    end subroutine 

    subroutine compute_mtrue_d3q19(rho, ux, uy, uz)
        use constants, only: zero 
        real(rkind), intent(in) :: rho, ux, uy, uz

        mvec_mrt_uniform(1 ) = rho
        mvec_mrt_uniform(2 ) = -11*rho + 19*rho*(ux*ux + uy*uy + uz*uz)
        mvec_mrt_uniform(3 ) = 3*rho - (11.d0/2.d0)*rho*(ux*ux + uy*uy + uz*uz)
        mvec_mrt_uniform(4 ) = rho*ux
        mvec_mrt_uniform(5 ) = -(2.d0/3.d0)*rho*ux
        mvec_mrt_uniform(6 ) = rho*uy
        mvec_mrt_uniform(7 ) = -(2.d0/3.d0)*rho*uy
        mvec_mrt_uniform(8 ) = rho*uz
        mvec_mrt_uniform(9 ) = -(2.d0/3.d0)*rho*uz
        mvec_mrt_uniform(10) = 2*rho*ux*ux - rho*uy*uy - rho*uz*uz
        mvec_mrt_uniform(11) = -rho*ux*ux + 0.5d0*uy*uy + 0.5d0*uz*uz
        mvec_mrt_uniform(12) = rho*uy*uy - rho*uz*uz
        mvec_mrt_uniform(13) = -0.5d0*rho*uy*uy + 0.5d0*rho*uz*uz
        mvec_mrt_uniform(14) = rho*ux*uy
        mvec_mrt_uniform(15) = rho*uz*uy
        mvec_mrt_uniform(16) = rho*ux*uz
        mvec_mrt_uniform(17) = zero
        mvec_mrt_uniform(18) = zero
        mvec_mrt_uniform(19) = zero
    end subroutine 

end module 

! These two subroutine are needed by d3q19 class, although they never get
! called because the problem is run in testing mode. 

subroutine initfields_bolt(decomp, inputfile, delta_x, rho, ux, uy, uz)
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info
    use constants, only: one, zero, two
    use d3q19_testing
    use exits, only: message
    use gridtools, only: linspace
    implicit none

    type(decomp_info), intent(in) :: decomp
    character(len=clen), intent(in) :: inputfile
    real(rkind), intent(in) :: delta_x
    real(rkind), dimension(:,:,:), intent(out) :: rho, ux, uy, uz  
    integer :: nz, idx

    select case(testing_type)
    case(0) ! Doesn't matter (just streaming test)
        rho = one
        ux = zero
        uy = zero
        uz = zero

        call message(0, "Velocity initialized for TEST 0 (streaming)")
    case(1)! Uniform flow (non-periodic)
        allocate(utrue(size(ux,1),size(ux,2),size(ux,3)))    
        allocate(vtrue(size(ux,1),size(ux,2),size(ux,3)))    
        allocate(wtrue(size(ux,1),size(ux,2),size(ux,3)))    
        rho = one
        ux = one
        uy = -one
        uz = one

        utrue = one
        vtrue = -one
        wtrue = one
        
        call compute_mtrue_d3q19(rho(1,1,1), ux(1,1,1), uy(1,1,1), uz(1,1,1))
        call message(0, "Velocity initialized for TEST 1 (uniform flow)")

    case(2) ! Laminar channel 
        allocate(utrue(size(ux,1),size(ux,2),size(ux,3)))    
        allocate(vtrue(size(ux,1),size(ux,2),size(ux,3)))    
        allocate(wtrue(size(ux,1),size(ux,2),size(ux,3)))    
       
        nz = size(ux,3)

        rho = one
        ux = one
        uy = zero
        uz = zero

        allocate(z(nz))
        z = linspace(-one,one,nz)
        do idx = 1,nz
            utrue(:,:,idx) = (Force*Reynolds_number/two)*(one - z(idx)*z(idx))
         end do 
        vtrue = zero
        wtrue = zero
        
        deallocate(z)
        call message(0, "Velocity initialized for TEST 2 (laminar channel)")
    case(3)
        allocate(utrue(size(ux,1),size(ux,2),size(ux,3)))    
        allocate(vtrue(size(ux,1),size(ux,2),size(ux,3)))    
        allocate(wtrue(size(ux,1),size(ux,2),size(ux,3)))    

        nz = size(ux,3)
        
        rho = one
        ux = one
        uy = zero
        uz = zero

        allocate(z(nz))
        z = linspace(-one,one,nz)
       
        call update_utrue_wormersley(0.d0)
        vtrue = zero
        wtrue = zero
      
        ux = utrue

        call message(0, "Velocity initialized for TEST 3 (Womersley flow)")
    end select 

end subroutine 

subroutine getWallBC_bolt(decomp, Re, delta_u, ux, uy, uz, uxB, uyB, uzB, uxT, uyT, uzT)
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info
    use constants, only: zero
    use d3q19_testing
    implicit none

    type(decomp_info), intent(in) :: decomp
    real(rkind), intent(in) :: Re, delta_u
    real(rkind), dimension(:,:,:), intent(in) :: ux, uy, uz
    real(rkind), dimension(:,:), intent(out) :: uxB, uyB, uzB, uxT, uyT, uzT
    integer :: nz

    nz = size(ux,3)

    uxB = utrue(:,:,1)
    uyB = vtrue(:,:,1) 
    uzB = wtrue(:,:,1) 

    uxT = utrue(:,:,nz)
    uyT = vtrue(:,:,nz)
    uzT = wtrue(:,:,nz)

end subroutine 


subroutine getWall_nut(decomp,delta_nu, ux, uy, uz, Re, tau_B, tau_T)
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info
    use constants, only: zero
    use d3q19_testing
    implicit none

    type(decomp_info), intent(in) :: decomp
    real(rkind), intent(in) :: delta_nu
    real(rkind), dimension(:,:,:), intent(in) :: ux, uy, uz
    real(rkind), intent(in) :: Re
    real(rkind), dimension(:,:), intent(out) :: tau_B, tau_T 
    integer :: nz

    nz = size(ux,3)

    tau_B = zero
    tau_T = zero

end subroutine 

subroutine getBodyForce(decomp, time, delta_u, delta_t, ux, uy, uz, Fx, Fy, Fz)
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info
    
    implicit none
    type(decomp_info), intent(in) :: decomp
    real(rkind), intent(in) :: time, delta_t, delta_u 
    real(rkind), dimension(:,:,:), intent(in) :: ux, uy, uz
    real(rkind), dimension(:,:,:), intent(out) :: Fx, Fy, Fz 


    Fx = 0.d0 
    Fy = 0.d0
    Fz = 0.d0 
end subroutine 



