module d3q19_channel3D
    use kind_parameters, only: rkind
    use constants, only: one, pi, two
    implicit none

    real(rkind), parameter :: Lxwant = two*pi
    real(rkind), parameter :: Lywant = two*pi
    real(rkind) :: utau  = 1.d0, uMatch 
    real(rkind) :: zmatch, zfirst
    real(rkind), dimension(:), allocatable :: x,y,z
    integer :: nx, ny, nz
    real(rkind) :: Lx, Ly, mfact

    ! Body force constants
    real(rkind), parameter :: kx_force = 1.d0
    real(rkind), parameter :: lambda1_force = 5.d0 
    real(rkind), parameter :: lambda2_force = 10.d0 
    real(rkind), parameter :: Fbase = 1.7150008761831d0
    real(rkind), parameter :: Force_amp = 15.d0     
end module 

module get_initial_profiles_channel
    use kind_parameters, only:rkind
    use constants, only: pi
    
    implicit none
    private
    public :: get_prof, get_musker_profile

    real(rkind), parameter :: Retau = 5.185915054711403d3
    real(rkind), parameter :: utau = 1.309580419899106d0  

contains

    pure subroutine get_prof(x,y,z,Lx,Ly,ux,uy,uz)
        real(rkind), intent(in) :: x,y,z, Lx, Ly
        real(rkind), intent(out) :: ux, uy, uz
        real(rkind) :: um
        integer :: nmodes_x, nmodes_y
        
        nmodes_x = 4
        nmodes_y = 2

        um = utau*get_musker_profile((1.d0-abs(z))*Retau)
        ux = um + get_womersley(z)*cos(2.d0*pi*nmodes_x*y/Ly)
        uy = get_womersley(z)*sin(2.d0*pi*nmodes_y*x/Lx)
        uz = 0.d0 

        !ux = um
        !uy = 0.d0
        
    end subroutine 

    pure elemental function get_womersley(z) result(u)
        use constants, only: imi
        real(rkind), intent(in) :: z
        real(rkind) :: u
        real(rkind) :: nu,omega,t, p0, D, alpha

        D = 1.d0; nu = 1.d-3; omega = 0.15d0; t=0.5d0; p0 = 1.d0;
        alpha = (D/2.d0)*sqrt(omega/nu)

        u = real((p0/(imi*omega))*(1.d0 - cosh(sqrt(2.d0)*(1 + imi)*alpha*z/D)/cosh(sqrt(2.d0)*(1 + &
          & imi)*alpha))*exp(imi*omega*t),rkind)
        
    end function 

    pure elemental function get_musker_profile(yp) result(up)
        real(rkind), intent(in) :: yp
        real(rkind) :: up

        up = 5.424d0*atan(0.11976047904*yp - 0.48802395209580833) + 0.434d0*log(((yp + 10.6)**9.6)/((yp**2 &
            & - 8.15*yp + 86)**2)) - 3.50727901936264842

    end function


end module 

module wall_module_routines
    use kind_parameters, only: rkind

contains
    pure function find_utau(utau_guess, umatch, zmatch, Re) result(utau)
        use get_initial_profiles_channel, only: get_musker_profile
        real(rkind), intent(in) :: utau_guess, umatch, zmatch, Re
        real(rkind) :: zplus, utau_diff, up1, utaunew
        real(rkind) :: utau

        utau_diff = utau_guess
        utau = utau_guess
        do while (utau_diff > 1.d-8)
            zplus = zmatch*utau*Re 
            up1 = get_musker_profile(zplus)                
            utaunew = umatch/up1
            utau_diff = abs(utaunew - utau)
            utau = utaunew
        end do 

        ! print*, tau
    end function 

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
    use exits, only: GracefulExit
    use get_initial_profiles_channel, only: get_prof

    implicit none

    type(decomp_info), intent(in) :: decomp
    character(len=clen), intent(in) :: inputfile
    real(rkind), intent(in) :: delta_x
    real(rkind), dimension(:,:,:), intent(out) :: rho, ux, uy, uz  
    integer :: idx, i, j, k, ii, jj 
    real(rkind), dimension(:), allocatable :: zE

    nx = decomp%xsz(1)
    ny = decomp%ysz(2)
    nz = decomp%zsz(3)
    
    allocate(zE(nz+1))
    zE = linspace(-1.d0,1.d0,nz+1)
    allocate(z(nz))
    z = 0.5d0*(zE(2:nz+1) + zE(1:nz))

    allocate(x(nx), y(ny))

    Lx = ceiling(Lxwant/delta_x)*delta_x
    Ly = ceiling(Lywant/delta_x)*delta_x

    if ((ceiling((Lxwant)/delta_x)).ne. nx) then
        print*, (ceiling((Lxwant)/delta_x)-1)
        print*, "Lx", Lx
        print*, "delta_x", delta_x
        print*, "nx", nx
        call GracefulExit("Inconsistent nx, ny, nz or delta_x provided",14)
    end if 

    if ((ceiling((Lywant)/delta_x)).ne. ny) then
        print*, (ceiling((Lywant)/delta_x)-1)
        print*, "Ly"
        call GracefulExit("Inconsistent nx, ny, nz or delta_x provided",14)
    end if 
   
    x = linspace(0.d0,Lx-delta_x,nx)
    y = linspace(0.d0,Ly-delta_x,ny)
    
    rho = one

    do k = decomp%zst(3),decomp%zen(3)
        jj = 1
        do j = decomp%zst(2),decomp%zen(2)
            ii = 1
            do i = decomp%zst(1),decomp%zen(1)
                call get_prof(x(i),y(j),z(k),Lx,Ly,ux(ii,jj,k),uy(ii,jj,k),uz(ii,jj,k))
                ii = ii + 1
            end do 
            jj = jj + 1
        end do 
    end do

    ux = 1.d0
    uy = 0.d0
    uz = 0.d0 

    zmatch = z(2) + 1.d0 
    zfirst = z(1) + 1.d0 
    mfact = 1.d0/real(nx*ny,rkind)

end subroutine 

subroutine getWallBC_bolt(decomp, Re, delta_u, ux, uy, uz, uxB, uyB, uzB, uxT, uyT, uzT)
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info, nrank 
    use constants, only: zero
    use d3q19_channel3D
    use reductions, only: p_sum 
    use wall_module_routines 
    use get_initial_profiles_channel, only: get_musker_profile
    use constants, only: zero
    use mpi 
    implicit none

    type(decomp_info), intent(in) :: decomp
    real(rkind), intent(in) :: Re, delta_u
    real(rkind), dimension(:,:,:), intent(in) :: ux, uy, uz
    real(rkind), dimension(:,:), intent(out) :: uxB, uyB, uzB, uxT, uyT, uzT
    real(rkind) :: utau_new, ubc
    integer :: ierr

    ! STEP 1: Get utau
    umatch = 0.5d0*mfact*(p_sum(sum(ux(:,:,2))) + p_sum(sum(ux(:,:,nz-1))))*delta_u
    utau_new = find_utau(utau, umatch, zmatch, Re) 
    utau = utau_new 
   
    ! STEP 2: Get ubc
    ubc = utau*get_musker_profile(zfirst*utau*Re)                

    ! STEP 3: Compute BC for bottom
    uxB = ux(:,:,2)
    uyB = uy(:,:,2)
    uzB = sqrt(uxB*uxB + uyB*uyB)
    uxB = uxB*ubc/(uzB*delta_u + 1.d-18)
    uyB = uyB*ubc/(uzB*delta_u + 1.d-18)
    uzB = zero

    ! STEP 3: Compute BC for top
    uxT = ux(:,:,nz-1)
    uyT = uy(:,:,nz-1)
    uzT = sqrt(uxT*uxT + uyT*uyT)
    uxT = uxT*ubc/(uzT*delta_u + 1.d-18)
    uyT = uyT*ubc/(uzT*delta_u + 1.d-18)
    uzT = zero
    
    !uxB = 0.d0
    !uyB = 0.d0
    !uzB = 0.d0

    !uxT = 0.d0
    !uyT = 0.d0
    !uzT = 0.d0

end subroutine 


subroutine getWall_nut(decomp, delta_nu, ux, uy, uz, Re, tau_B, tau_T)
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info
    use constants, only: zero
    use d3q19_channel3D
    use constants, only: three, half
    use exits, only: message
    implicit none

    type(decomp_info), intent(in) :: decomp
    real(rkind), intent(in) :: delta_nu
    real(rkind), dimension(:,:,:), intent(in) :: ux, uy, uz
    real(rkind), intent(in) :: Re
    real(rkind), dimension(:,:), intent(out) :: tau_B, tau_T 
    real(rkind) :: yp, dudy_wall, Va, nu_t
    real(rkind), parameter :: kappa = 0.384d0

    yp = zfirst*utau*Re
    dudy_wall = (2939515020475289880216614001228638060544.d0/(125.d0*(519229685836867467736034237096533489.d0*yp*yp &
              & - 4231721939638177189649975113492201472.d0*yp + 44824125224070191636923425342095360000.d0)) + &
              & (60760.d0*yp*yp - 1132089.d0*yp + 10832423.d0)/(250.d0*(245.d0*yp*yp - 39.d0*yp + 100.d0*yp*yp*yp + 91160.d0)))

    dudy_wall = dudy_wall*utau*utau*Re
    Va = 1.d0 - exp(-yp/26.d0)
    nu_t = (((kappa*zfirst*Va)**2)*dudy_wall + (one/Re))/delta_nu

    call message(2,"Wall nu_t:", nu_t)
    tau_B = nu_t*three + half 
    tau_T = nu_t*three + half 

end subroutine 

subroutine getBodyForce(decomp, time, delta_u, delta_t, ux, uy, uz, Fx, Fy, Fz)
    use d3q19_channel3D
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info
    use constants, only: pi
    
    implicit none
    type(decomp_info), intent(in) :: decomp
    real(rkind), intent(in) :: time, delta_t, delta_u 
    real(rkind), dimension(:,:,:), intent(in) :: ux, uy, uz
    real(rkind), dimension(:,:,:), intent(out) :: Fx, Fy, Fz 
    real(rkind) :: Fconst, fact_time, zfact
    integer :: i, j, k

    Fconst = Fbase*(1.d0 - exp(-lambda2_force*time))
    fact_time = exp(-lambda1_force*time)

    do k = decomp%zst(3),decomp%zen(3)
        zfact = Force_amp*(tan(z(k)) + 0.2d0*cos(1.d0*z(k)*2.d0*pi))
        do j = decomp%zst(2),decomp%zen(2)
            !$omp simd
            do i = decomp%zst(1),decomp%zen(1)
                Fx(i-decomp%zst(1)+1,j-decomp%zst(2)+1,k) = (Fconst +  zfact*fact_time*sin(kx_force*x(i)*2.d0*pi/Lx))*delta_t/delta_u
                Fy(i-decomp%zst(1)+1,j-decomp%zst(2)+1,k) = (zfact*fact_time*sin(kx_force*y(j)*2.d0*pi/Ly))*delta_t/delta_u
            end do 
        end do 
    end do 
    Fz = 0.d0 
end subroutine 



