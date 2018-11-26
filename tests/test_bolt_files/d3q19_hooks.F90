module d3q19_testing
    use kind_parameters, only: rkind
    implicit none
    integer :: testing_type = 0
    real(rkind) :: Reynolds_number = 20._rkind, Force = 1._rkind
    real(rkind), dimension(:,:,:), allocatable :: utrue, vtrue, wtrue
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
    real(rkind), dimension(:), allocatable :: z
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
        uy = one
        uz = one

        utrue = one
        vtrue = one
        wtrue = one
        
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
            utrue(:,:,idx) = (Reynolds_number/two)*(one - z(idx)*z(idx))
         end do 
        vtrue = zero
        wtrue = zero
        
        deallocate(z)

        call message(0, "Velocity initialized for TEST 2 (laminar channel)")
    end select 

end subroutine 

subroutine getWallBC_bolt(decomp, ux, uy, uz, uxB, uyB, uzB, uxT, uyT, uzT)
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info
    use constants, only: zero
    use d3q19_testing
    implicit none

    type(decomp_info), intent(in) :: decomp
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
