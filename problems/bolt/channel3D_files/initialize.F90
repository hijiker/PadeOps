module d3q19_channel3D
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
    real(rkind), parameter :: Force_amp = 15.d0     
    
    real(rkind), parameter :: umean_target = 0.1d0
    real(rkind), parameter :: uwall = 0.d0, vwall = 0.d0, wwall = 0.d0 

    real(rkind) :: umean = 0.d0  
    integer, parameter :: utau_method = 1
    real(rkind), dimension(:,:), allocatable :: utau_up, utau_do, yp_up, yp_do, dudz_up, dudz_do, ubc_up, ubc_do

    integer, parameter :: NoslipBC     = 0
    integer, parameter :: LogLawWMBC = 1
    integer, parameter :: FullLogLawWMBC = 2
    integer, parameter :: Forcing_Type = 1

    integer            :: BC_type = LogLawWMBC !NoslipBC
    real(rkind)        :: Fbase = 2.953879d-03 

contains
    subroutine allocate_WM_arrays(decomp)
        use decomp_2d, only: decomp_info
        type(decomp_info), intent(in) :: decomp

        allocate(utau_up(decomp%zsz(1),decomp%zsz(2)))
        allocate(utau_do(decomp%zsz(1),decomp%zsz(2)))
        allocate(yp_up(decomp%zsz(1),decomp%zsz(2)))
        allocate(yp_do(decomp%zsz(1),decomp%zsz(2)))
        
        allocate(dudz_up(decomp%zsz(1),decomp%zsz(2)))
        allocate(dudz_do(decomp%zsz(1),decomp%zsz(2)))
        
        allocate(ubc_up(decomp%zsz(1),decomp%zsz(2)))
        allocate(ubc_do(decomp%zsz(1),decomp%zsz(2)))
        
        utau_up = 1.d0 
        utau_do = 1.d0 
    end subroutine 
end module 

module get_initial_profiles_channel
    use kind_parameters, only:rkind
    use constants, only: pi
    
    implicit none
    private
    public :: Retau, utau, get_prof, get_musker_profile, get_musker_gradient

    real(rkind) :: Retau = 5.434960d02 
    real(rkind) :: utau = 5.434960d-02 


contains

    pure subroutine get_prof(xG,yG,zG,ux,uy,uz)
        real(rkind), intent(in) :: xG,yG,zG
        real(rkind), intent(out) :: ux, uy, uz
        real(rkind) :: Lx, Ly, um, nmodes, lambda, Amp 

        um = utau*get_musker_profile((1.d0-abs(zG))*Retau)
        nmodes = 2
        lambda = 0.5d0
        Amp = 0.05d0
        Lx = 6.d0
        Ly = 3.d0 

        uz = Amp*((nmodes*2*pi/Lx)*exp(-lambda*(zG - 1)*(zG + 1))*((zG - 1)**2)*((zG + 1)**2)*sin(nmodes*xG*2*pi/Lx))*sin(nmodes*yG*2*pi/Ly)
        ux = um + Amp*(2*zG*exp(-lambda*(zG**2 - 1))*cos(nmodes*xG*2*pi/Lx)*sin(nmodes*yG*2*pi/Ly)*(zG**2 - 1)*(lambda - lambda*zG**2 + 2))
        uy = 0.d0
        
    end subroutine 

    pure elemental function get_musker_profile(yp) result(up)
        real(rkind), intent(in) :: yp
        real(rkind) :: up

        up = 5.424d0*atan(0.11976047904*yp - 0.48802395209580833) + 0.434d0*log(((yp + 10.6)**9.6)/((yp**2 &
            & - 8.15*yp + 86)**2)) - 3.50727901936264842

    end function
    
    pure elemental function get_musker_gradient(yp) result(dudy_wall)
        real(rkind), intent(in) :: yp
        real(rkind) :: dudy_wall

        dudy_wall = (2939515020475289880216614001228638060544.d0/(125.d0*(519229685836867467736034237096533489.d0*yp*yp &
              & - 4231721939638177189649975113492201472.d0*yp + 44824125224070191636923425342095360000.d0)) + &
              & (60760.d0*yp*yp - 1132089.d0*yp + 10832423.d0)/(250.d0*(245.d0*yp*yp - 39.d0*yp + 100.d0*yp*yp*yp + 91160.d0)))

    end function


end module 

module wall_model_routines
    use kind_parameters, only: rkind
    implicit none
contains
    pure function find_utau(utau_guess, umatch, zmatch, Re) result(ustar)
        real(rkind), intent(in) :: utau_guess, umatch, zmatch, Re
        real(rkind) :: ustar, g, dgdus
        integer, parameter :: iter_lim = 50
        integer :: iter

        ustar = utau_guess
        iter = 0
        call g_ustar(ustar,umatch,zmatch,Re,g,dgdus)
        do while ((abs(g) > 1.d-8) .and. (iter<iter_lim))
            ustar = ustar - g/(dgdus+1E-14)
            call g_ustar(ustar,umatch,zmatch,Re,g,dgdus)
            iter = iter + 1
        end do 
    end function 

    pure subroutine g_ustar(ustar,umatch,zmatch,Re,g,dgdus)
        use get_initial_profiles_channel, only: get_musker_profile, get_musker_gradient
        real(rkind), intent(in) :: ustar, umatch, zmatch, Re
        real(rkind), intent(out) :: g, dgdus
        real(rkind) :: yp, up, dupdyp

        yp = zmatch*Re*ustar
        up = get_musker_profile(yp)
        dupdyp = get_musker_gradient(yp)
        
        g = umatch - ustar*up
        dgdus = -up - ustar*Re*zmatch*dupdyp
    end subroutine 
end module 


! These two subroutine are needed by d3q19 class, although they never get
! called because the problem is run in testing mode. 
subroutine initfields_bolt(decomp, inputfile, delta_x, rho, ux, uy, uz)
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info, nrank 
    use constants, only: one, zero, two
    use d3q19_channel3D
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
    real(rkind), dimension(:), allocatable :: zE
    real(rkind), parameter :: Noise_Amp = 1.d2
    real(rkind), dimension(:,:,:), allocatable :: randArr
    integer :: seedu = 1394, seedv = 6322, seedw = 7543

    nx = decomp%xsz(1)
    ny = decomp%ysz(2)
    nz = decomp%zsz(3)
    
    allocate(zE(nz+1))
    zE = linspace(-1.d0,1.d0,nz+1)
    allocate(z(nz))
    z = 0.5d0*(zE(2:nz+1) + zE(1:nz))

    allocate(x(nx), y(ny))

    Lx = (real(nx,rkind)/real(nz,rkind))*2.d0 
    Ly = (real(ny,rkind)/real(nz,rkind))*2.d0 
   
    x = linspace(0.d0,Lx-delta_x,nx)
    y = linspace(0.d0,Ly-delta_x,ny)
    
    if (abs((x(2) - x(1)) - (z(2) - z(1))) > 1.d-14) then
        print*, "x_issue", (x(2) - x(1)), (z(2) - z(1)), delta_x
        call gracefulExit("Incorrect choice of nx, nz and delta_x",13)
    end if 
    if (abs((y(2) - y(1)) - (z(2) - z(1))) > 1.d-14) then
        print*, "y_issue", (y(2) - y(1)), (z(2) - z(1)), delta_x
        call gracefulExit("Incorrect choice of ny, nz and delta_x",14)
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

    zmatch = z(2) + 1.d0 
    zfirst = z(1) + 1.d0 
    mfact = 1.d0/real(nx*ny,rkind)
        
    call allocate_WM_arrays(decomp)

end subroutine 


subroutine getWallBC_bolt(decomp, Re, delta_u, ux, uy, uz, uxB, uyB, uzB, uxT, uyT, uzT)
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info
    use constants, only: zero
    use d3q19_channel3D
    use reductions, only: p_sum, p_maxval
    use wall_model_routines 
    use get_initial_profiles_channel, only: Retau, utau, get_musker_profile
    use constants, only: zero
    use mpi 
    
    implicit none

    type(decomp_info), intent(in) :: decomp
    real(rkind), intent(in) :: Re, delta_u
    real(rkind), dimension(:,:,:), intent(in) :: ux, uy, uz
    real(rkind), dimension(:,:), intent(out) :: uxB, uyB, uzB, uxT, uyT, uzT
    real(rkind) :: utau_new, umatchB, umatchT, onebydelta, uplus_match, ex2, ey2
    integer :: i, j

    nz = decomp%zsz(3)
    select case (BC_Type)
    case (FullLogLawWMBC)
        if (.not. allocated(utau_up)) then
            call allocate_WM_arrays(decomp)
        end if 
        
        onebydelta = one/delta_u
        
        do j = 1,decomp%zsz(2)
            do i = 1,decomp%zsz(1)
                
                umatchB = delta_u*sqrt(ux(i,j,2)*ux(i,j,2) + uy(i,j,2)*uy(i,j,2))
                utau_new = find_utau(utau_do(i,j), umatchB, zmatch, Re) 
                utau_do(i,j) = utau_new 
                
                umatchT = delta_u*sqrt(ux(i,j,nz-1)*ux(i,j,nz-1) + uy(i,j,nz-1)*uy(i,j,nz-1))
                utau_new = find_utau(utau_up(i,j), umatchT, zmatch, Re) 
                utau_up(i,j) = utau_new 
                
                ubc_do(i,j) = utau_do(i,j)*get_musker_profile(zfirst*utau_do(i,j)*Re)                
                ubc_up(i,j) = utau_up(i,j)*get_musker_profile(zfirst*utau_up(i,j)*Re)           
                
                ex2 = ux(i,j,2)*delta_u/umatchB
                ey2 = uy(i,j,2)*delta_u/umatchB
                uxB(i,j) = (onebydelta*ubc_do(i,j))*ex2
                uyB(i,j) = (onebydelta*ubc_do(i,j))*ey2
                
                ex2 = ux(i,j,nz-1)*delta_u/umatchT
                ey2 = uy(i,j,nz-1)*delta_u/umatchT
                uxT(i,j) = (onebydelta*ubc_up(i,j))*ex2
                uyT(i,j) = (onebydelta*ubc_up(i,j))*ey2

            end do 
        end do 

        uzB = (8.d0/15.d0)*(wwall/delta_u + (5.d0/4.d0)*uz(:,:,2) - (3.d0/8.d0)*uz(:,:,3))
        uzT = (8.d0/15.d0)*(wwall/delta_u + (5.d0/4.d0)*uz(:,:,nz-1) - (3.d0/8.d0)*uz(:,:,nz-2))
    
    case (LogLawWMBC)
        if (.not. allocated(utau_up)) then
            call allocate_WM_arrays(decomp)
        end if 

        onebydelta = one/delta_u
        
        do j = 1,decomp%zsz(2)
            do i = 1,decomp%zsz(1)
                
                uplus_match = get_musker_profile(zmatch*Retau) + 1.d-14
                
                umatchB = delta_u*sqrt(ux(i,j,2)*ux(i,j,2) + uy(i,j,2)*uy(i,j,2))
                utau_do(i,j) = umatchB/uplus_match
                
                umatchT = delta_u*sqrt(ux(i,j,nz-1)*ux(i,j,nz-1) + uy(i,j,nz-1)*uy(i,j,nz-1))
                utau_up(i,j) = umatchT/uplus_match
                
                ubc_do(i,j) = utau_do(i,j)*get_musker_profile(zfirst*Retau)                
                ubc_up(i,j) = utau_up(i,j)*get_musker_profile(zfirst*Retau)           
                
                ex2 = ux(i,j,2)*delta_u/umatchB
                ey2 = uy(i,j,2)*delta_u/umatchB
                uxB(i,j) = (onebydelta*ubc_do(i,j))*ex2
                uyB(i,j) = (onebydelta*ubc_do(i,j))*ey2
                
                ex2 = ux(i,j,nz-1)*delta_u/umatchT
                ey2 = uy(i,j,nz-1)*delta_u/umatchT
                uxT(i,j) = (onebydelta*ubc_up(i,j))*ex2
                uyT(i,j) = (onebydelta*ubc_up(i,j))*ey2

            end do 
        end do 
        

        uzB = (8.d0/15.d0)*(wwall/delta_u + (5.d0/4.d0)*uz(:,:,2) - (3.d0/8.d0)*uz(:,:,3))
        uzT = (8.d0/15.d0)*(wwall/delta_u + (5.d0/4.d0)*uz(:,:,nz-1) - (3.d0/8.d0)*uz(:,:,nz-2))

    case (NoSlipBC)
        uxB = (8.d0/15.d0)*(uwall/delta_u + (5.d0/4.d0)*ux(:,:,2) - (3.d0/8.d0)*ux(:,:,3))
        uyB = (8.d0/15.d0)*(vwall/delta_u + (5.d0/4.d0)*uy(:,:,2) - (3.d0/8.d0)*uy(:,:,3))
        uzB = (8.d0/15.d0)*(wwall/delta_u + (5.d0/4.d0)*uz(:,:,2) - (3.d0/8.d0)*uz(:,:,3))
        
        uxT = (8.d0/15.d0)*(uwall/delta_u + (5.d0/4.d0)*ux(:,:,nz-1) - (3.d0/8.d0)*ux(:,:,nz-2))
        uyT = (8.d0/15.d0)*(vwall/delta_u + (5.d0/4.d0)*uy(:,:,nz-1) - (3.d0/8.d0)*uy(:,:,nz-2))
        uzT = (8.d0/15.d0)*(wwall/delta_u + (5.d0/4.d0)*uz(:,:,nz-1) - (3.d0/8.d0)*uz(:,:,nz-2))
    end select 

end subroutine 


subroutine getWall_nut(decomp, delta_nu, ux, uy, uz, Re, tau_B, tau_T)
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info
    use constants, only: zero
    use d3q19_channel3D
    use constants, only: three, half
    use exits, only: message
    use get_initial_profiles_channel, only: get_musker_gradient
    implicit none

    type(decomp_info), intent(in) :: decomp
    real(rkind), intent(in) :: delta_nu
    real(rkind), dimension(:,:,:), intent(in) :: ux, uy, uz
    real(rkind), intent(in) :: Re
    real(rkind), dimension(:,:), intent(out) :: tau_B, tau_T 
    real(rkind) ::nu, oneByCsq
    
    oneByCsq = 3.d0 
    nu = (one/Re)/delta_nu
    tau_B = half + oneByCsq*nu
    tau_T = half + oneByCsq*nu

   

end subroutine 

subroutine getBodyForce(decomp, time, delta_u, delta_t, ux, uy, uz, Fx, Fy, Fz)
    use d3q19_channel3D
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info
    use constants, only: pi
    use reductions, only: p_sum
    
    implicit none
    type(decomp_info), intent(in) :: decomp
    real(rkind), intent(in) :: time, delta_t, delta_u 
    real(rkind), dimension(:,:,:), intent(in) :: ux, uy, uz
    real(rkind), dimension(:,:,:), intent(out) :: Fx, Fy, Fz 
    real(rkind) :: Fconst, fact_time, zfact
    integer :: i, j, k

    real(rkind) :: mnorm

    mnorm = 1.d0/(real(decomp%xsz(1),rkind)*real(decomp%ysz(2),rkind)*real(decomp%zsz(3),rkind))
    umean= p_sum(sum(ux))*mnorm

    select case (Forcing_type)
    case (1)
        Fz = 0.d0 
        Fy = 0.d0
        Fx = Fbase*delta_t/delta_u + delta_t*(umean_target - umean)
    case default
        Fz = 0.d0 
        Fy = 0.d0
        Fx = Fbase*delta_t/delta_u
    end select 

end subroutine 



