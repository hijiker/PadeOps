module d3q19_channel3D
    use kind_parameters, only: rkind
    use constants, only: one, pi, two
    implicit none

    real(rkind), parameter :: Lxwant = two*pi
    real(rkind), parameter :: Lywant = two*pi
    real(rkind) :: utau  
    real(rkind), dimension(:), allocatable :: x,y,z
    integer :: nx, ny, nz
    real(rkind) :: Lx, Ly, delta_x

end module 

! These two subroutine are needed by d3q19 class, although they never get
! called because the problem is run in testing mode. 

subroutine initfields_bolt(decomp, inputfile, delta_x, rho, ux, uy, uz)
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info
    use constants, only: one, zero, two
    use d3q19_channel3D
    use exits, only: message
    use gridtools, only: linspace
    implicit none

    type(decomp_info), intent(in) :: decomp
    character(len=clen), intent(in) :: inputfile
    real(rkind), intent(in) :: delta_x
    real(rkind), dimension(:,:,:), intent(out) :: rho, ux, uy, uz  
    integer :: nz, idx
    real(rkind), dimension(:), allocatable :: zE

    nx = decom%xsz(1)
    ny = decom%ysz(2)
    nz = decom%zsz(3)
    
    allocate(zE(nz+1))
    zE = linspace(-1.d0,1.d0,nz+1)
    allocate(z(nz))
    z = 0.5d0*(zE(2:nz+1) + zE(1:nz))

    allocate(x(nx), y(ny))
    delta_x = z(2) - z(1)

    Lx = ceiling(Lxwant/delta_x)*delta_x
    Lx = ceiling(Lywant/delta_x)*delta_x

    if ((ceiling((Lxwant)/delta_x)-1).ne. nx) then
        call GracefulExit("Inconsistent nx, ny, nz provided",14)
    end if 

    if ((ceiling((Lywant)/delta_x)-1).ne. ny) then
        call GracefulExit("Inconsistent nx, ny, nz provided",14)
    end if 
   
    x = linspace(0.d0,2.d0*pi-delta_x,nx)
    y = linspace(0.d0,2.d0*pi-delta_x,ny)
    
    rho = one
    ux = one
    uy = zero
    uz = zero

    allocate(z(nz))
    z = linspace(-one,one,nz)
    rho = one
    ux = one
    uy = zero
    uz = zero

    allocate(z(nz))
    z = linspace(-one,one,nz)
    x = linspace
    
    call update_utrue_wormersley(0.d0)
    vtrue = zero
    wtrue = zero
    
    ux = utrue

    call message(0, "Velocity initialized for Channel3D")

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


subroutine getWall_nut(decomp,ux, uy, uz, Re, tau_B, tau_T)
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info
    use constants, only: zero
    use d3q19_testing
    implicit none

    type(decomp_info), intent(in) :: decomp
    real(rkind), dimension(:,:,:), intent(in) :: ux, uy, uz
    real(rkind), intent(in) :: Re
    real(rkind), dimension(:,:), intent(out) :: tau_B, tau_T 
    integer :: nz

    nz = size(ux,3)

    tau_B = zero
    tau_T = zero

end subroutine 

