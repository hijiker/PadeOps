module d3q19_channelDNS
    use kind_parameters, only: rkind
    use constants, only: one, pi, two
    implicit none

    real(rkind) :: zmatch, zfirst
    real(rkind), dimension(:), allocatable :: x,y,z
    integer :: nx, ny, nz
    real(rkind) :: Lx, Ly, mfact

    ! Body force constants
    real(rkind), parameter :: kx_force = 1.d0
    real(rkind), parameter :: lambda1_force = 5.d0 
    real(rkind), parameter :: lambda2_force = 10.d0 
    real(rkind), parameter :: Fbase = 2.953879d-03 
    real(rkind), parameter :: Force_amp = 15.d0     
   
end module 

module get_initial_profiles_channel
    use kind_parameters, only:rkind
    use constants, only: pi, zero, two
    
    implicit none
    private
    public :: utau, get_prof, get_musker_profile

    real(rkind), parameter :: Retau = 5.434960d02 
    real(rkind), parameter :: utau = 5.434960d-02 

contains
    pure subroutine get_prof(xG,yG,zG,ux,uy,uz)
        real(rkind), intent(in) :: xG,yG,zG
        real(rkind), intent(out) :: ux, uy, uz
        real(rkind) :: Lx, um, nmodes, lambda, Amp 

        um = utau*get_musker_profile((1.d0-abs(zG))*Retau)
        nmodes = 2
        lambda = 0.5d0
        Amp = 0.05d0
        Lx = 6.d0 

        uz = Amp*((nmodes*2*pi/Lx)*exp(-lambda*(zG - 1)*(zG + 1))*((zG - 1)**2)*((zG + 1)**2)*sin(nmodes*xG*2*pi/Lx))
        ux = um + Amp*(2*zG*exp(-lambda*(zG**2 - 1))*cos(nmodes*xG*2*pi/Lx)*(zG**2 - 1)*(lambda - lambda*zG**2 + 2))
        uy = 0.d0
        
    end subroutine 

    pure elemental function get_musker_profile(yp) result(up)
        real(rkind), intent(in) :: yp
        real(rkind) :: up

        up = 5.424d0*atan(0.11976047904*yp - 0.48802395209580833) + 0.434d0*log(((yp + 10.6)**9.6)/((yp**2 &
            & - 8.15*yp + 86)**2)) - 3.50727901936264842

    end function
    
end module 


! These two subroutine are needed by d3q19 class, although they never get
! called because the problem is run in testing mode. 
subroutine initfields_bolt(decomp, inputfile, delta_x, rho, ux, uy, uz)
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info, nrank 
    use constants, only: one, zero, two
    use d3q19_channelDNS
    use exits, only: message
    use gridtools, only: linspace
    use exits, only: GracefulExit
    use get_initial_profiles_channel, only: get_prof
    use random,             only: gaussian_random

    implicit none

    type(decomp_info), intent(in) :: decomp
    character(len=clen), intent(in) :: inputfile
    real(rkind), intent(in) :: delta_x
    real(rkind), dimension(:,:,:), intent(out) :: rho, ux, uy, uz  
    integer :: i, j, k, ii, jj 
    real(rkind), parameter :: Noise_Amp = 1.d-3
    real(rkind), dimension(:,:,:), allocatable :: randArr
    integer :: seedu = 1394, seedv = 6322, seedw = 7543

    nx = decomp%xsz(1)
    ny = decomp%ysz(2)
    nz = decomp%zsz(3)
    
    allocate(z(nz))
    z = linspace(-1.d0,1.d0,nz) 

    allocate(x(nx), y(ny))

    Lx = (real(nx,rkind)/real(nz-1,rkind))*2.d0 
    Ly = (real(ny,rkind)/real(nz-1,rkind))*2.d0 
   
    x = linspace(0.d0,Lx-delta_x,nx)
    y = linspace(0.d0,Ly-delta_x,ny)
    
    if (abs((x(2) - x(1)) - (z(2) - z(1))) > 1.d-14) then
        call gracefulExit("Incorrect choice of nx, nz and delta_x",13)
    end if 
    if (abs((y(2) - y(1)) - (z(2) - z(1))) > 1.d-14) then
        call gracefulExit("Incorrect choice of ny, nz and delta_x",13)
    end if 
    
    rho = one

    do k = decomp%zst(3),decomp%zen(3)
        jj = 1
        do j = decomp%zst(2),decomp%zen(2)
            ii = 1
            do i = decomp%zst(1),decomp%zen(1)
                call get_prof(x(i),y(j),z(k),ux(ii,jj,k),uy(ii,jj,k),uz(ii,jj,k))
                ii = ii + 1
            end do 
            jj = jj + 1
        end do 
    end do
   
    allocate(randArr(size(ux,1),size(ux,2),size(ux,3)))
    call gaussian_random(randArr,zero,one,seedu + 100*nrank)
    ux  = ux + Noise_Amp*randArr
    
    call gaussian_random(randArr,zero,one,seedv + 100*nrank)
    uy  = uy + Noise_Amp*randArr
    
    call gaussian_random(randArr,zero,one,seedw + 100*nrank)
    uz  = uz + Noise_Amp*randArr
    deallocate(randArr)

        

end subroutine 


subroutine getWallBC_bolt(decomp, Re, delta_u, ux, uy, uz, uxB, uyB, uzB, uxT, uyT, uzT)
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info
    use constants, only: zero
    use d3q19_channelDNS
    use reductions, only: p_sum, p_maxval
    use get_initial_profiles_channel, only: utau, get_musker_profile
    use constants, only: zero
    use mpi 
    
    implicit none

    type(decomp_info), intent(in) :: decomp
    real(rkind), intent(in) :: Re, delta_u
    real(rkind), dimension(:,:,:), intent(in) :: ux, uy, uz
    real(rkind), dimension(:,:), intent(out) :: uxB, uyB, uzB, uxT, uyT, uzT
    real(rkind) :: utau_new, umatch, onebydelta
    integer :: i, j


    uxB = zero 
    uyB = zero
    uzB = zero

    uxT = zero 
    uyT = zero 
    uzT = zero 


end subroutine 


subroutine getWall_nut(decomp, delta_nu, ux, uy, uz, Re, tau_B, tau_T)
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info
    use constants, only: zero
    use d3q19_channelDNS
    use constants, only: half 
    use exits, only: message
    implicit none

    type(decomp_info), intent(in) :: decomp
    real(rkind), intent(in) :: delta_nu
    real(rkind), dimension(:,:,:), intent(in) :: ux, uy, uz
    real(rkind), intent(in) :: Re
    real(rkind), dimension(:,:), intent(out) :: tau_B, tau_T 

    tau_B = half
    tau_T = half

end subroutine 

subroutine getBodyForce(decomp, time, delta_u, delta_t, ux, uy, uz, Fx, Fy, Fz)
    use d3q19_channelDNS
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info
    use constants, only: pi
    
    implicit none
    type(decomp_info), intent(in) :: decomp
    real(rkind), intent(in) :: time, delta_t, delta_u 
    real(rkind), dimension(:,:,:), intent(in) :: ux, uy, uz
    real(rkind), dimension(:,:,:), intent(out) :: Fx, Fy, Fz 

    Fz = 0.d0 
    Fy = 0.d0
    Fx = Fbase*delta_t/delta_u
end subroutine 



