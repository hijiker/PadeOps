module d3q19_taylorgreen3D
    use kind_parameters, only: rkind
    use constants, only: one, pi, two
    implicit none

    integer, parameter :: TG_mode = 1 


    real(rkind) :: u0 = 2.00450365d0
    real(rkind), dimension(:), allocatable :: x,y,z
    integer :: nx, ny, nz
    ! Body force constants

end module 

subroutine initfields_bolt(decomp, inputfile, delta_x, rho, ux, uy, uz)
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info
    use constants, only: one, zero, two
    use d3q19_taylorgreen3D
    use exits, only: message
    use gridtools, only: linspace
    use exits, only: GracefulExit
    implicit none

    type(decomp_info), intent(in) :: decomp
    character(len=clen), intent(in) :: inputfile
    real(rkind), intent(in) :: delta_x
    real(rkind), dimension(:,:,:), intent(out) :: rho, ux, uy, uz  
    integer :: idx, i, j, k, ii, jj 
    real(rkind) :: delta_true


    nx = decomp%xsz(1)
    ny = decomp%ysz(2)
    nz = decomp%zsz(3)

    allocate(x(nx), y(ny), z(nz))
 
    delta_true = two*pi/real(nx,rkind)
    if (abs(delta_x - delta_true)>1.D-14) then
        print*, "Delta_x", delta_x
        print*, "True", delta_true
        call gracefulExit("Incorrect delta_x prescribed",23)
    end if 
    x = linspace(0.d0,two*pi-delta_x,nx) - pi 
    y = linspace(0.d0,two*pi-delta_x,ny) - pi
    z = linspace(0.d0,two*pi-delta_x,nz) - pi
    
    rho = one

    select case(TG_mode)
    case(1)
        do k = decomp%zst(3),decomp%zen(3)
            jj = 1
            do j = decomp%zst(2),decomp%zen(2)
                ii = 1
                do i = decomp%zst(1),decomp%zen(1)
                    ux(ii,jj,k) =  u0*sin(x(i))*cos(y(j))*cos(z(k))  
                    uy(ii,jj,k) = -u0*cos(x(i))*sin(y(j))*cos(z(k))  
                    ii = ii + 1
                end do 
                jj = jj + 1
            end do 
        end do
        uz = 0.d0 
    case(2)
        do k = decomp%zst(3),decomp%zen(3)
            jj = 1
            do j = decomp%zst(2),decomp%zen(2)
                ii = 1
                do i = decomp%zst(1),decomp%zen(1)
                    ii = ii + 1
                end do 
                jj = jj + 1
            end do 
        end do
        uy = 0.d0 
    case(3)
        do k = decomp%zst(3),decomp%zen(3)
            jj = 1
            do j = decomp%zst(2),decomp%zen(2)
                ii = 1
                do i = decomp%zst(1),decomp%zen(1)
                    ii = ii + 1
                end do 
                jj = jj + 1
            end do 
        end do
        ux = 0.d0 
    end select 

end subroutine 

subroutine getWallBC_bolt(decomp, Re, delta_u, ux, uy, uz, uxB, uyB, uzB, uxT, uyT, uzT)
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info, nrank 
    use constants, only: zero
    use reductions, only: p_sum 
    use constants, only: zero
    use mpi 
    implicit none

    type(decomp_info), intent(in) :: decomp
    real(rkind), intent(in) :: Re, delta_u
    real(rkind), dimension(:,:,:), intent(in) :: ux, uy, uz
    real(rkind), dimension(:,:), intent(out) :: uxB, uyB, uzB, uxT, uyT, uzT

    ! nothing to do here
    uxB = 0.d0
    uyB = 0.d0
    uzB = 0.d0
    
    uxT = 0.d0
    uyT = 0.d0
    uzT = 0.d0

end subroutine 


subroutine getWall_nut(decomp, delta_nu, ux, uy, uz, Re, tau_B, tau_T)
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info
    use exits, only: message
    implicit none

    type(decomp_info), intent(in) :: decomp
    real(rkind), intent(in) :: delta_nu
    real(rkind), dimension(:,:,:), intent(in) :: ux, uy, uz
    real(rkind), intent(in) :: Re
    real(rkind), dimension(:,:), intent(out) :: tau_B, tau_T 

    tau_B = 0.d0
    tau_T = 0.d0 

end subroutine 

subroutine getBodyForce(decomp, time, delta_u, delta_t, ux, uy, uz, Fx, Fy, Fz)
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info
    use constants, only: pi
    
    implicit none
    type(decomp_info), intent(in) :: decomp
    real(rkind), intent(in) :: time, delta_t, delta_u 
    real(rkind), dimension(:,:,:), intent(in) :: ux, uy, uz
    real(rkind), dimension(:,:,:), intent(out) :: Fx, Fy, Fz 

    Fx = 0.d0
    Fy = 0.d0
    Fz = 0.d0 

end subroutine 



