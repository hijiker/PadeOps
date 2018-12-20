module d3q19_taylorgreen2D
    use kind_parameters, only: rkind
    use constants, only: one, pi, two
    implicit none

    integer, parameter :: TG_mode = 3 


    real(rkind), parameter :: Lx = two*pi
    real(rkind), parameter :: Ly = two*pi
    real(rkind), parameter :: Lz = two*pi
    real(rkind) :: utau  = 1.d0, uMatch 
    real(rkind), dimension(:), allocatable :: x,y,z
    integer :: nx, ny, nz
    real(rkind), dimension(:,:,:), allocatable :: uexact, vexact, wexact
    ! Body force constants

contains

    subroutine update_exact_TG(time, Re, decomp)
    use decomp_2d, only: decomp_info
        real(rkind), intent(in) :: time, Re
        integer :: i, j, k, ii, jj, kk
        type(decomp_info), intent(in) :: decomp
        real(rkind) :: Ft

        if (.not. allocated(uexact)) allocate(uexact(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3))) 
        if (.not. allocated(vexact)) allocate(vexact(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3))) 
        if (.not. allocated(wexact)) allocate(wexact(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3))) 
        
        Ft = exp(-2.d0*time/Re)

        select case(TG_mode)
        case(1)
            do k = decomp%zst(3),decomp%zen(3)
                jj = 1
                do j = decomp%zst(2),decomp%zen(2)
                    ii = 1
                    do i = decomp%zst(1),decomp%zen(1)
                        uexact(ii,jj,k) =  cos(x(i))*sin(y(j))*Ft
                        vexact(ii,jj,k) = -sin(x(i))*cos(y(j))*Ft
                        ii = ii + 1
                    end do 
                    jj = jj + 1
                end do 
            end do
            wexact = 0.d0 
        case(2)
            do k = decomp%zst(3),decomp%zen(3)
                jj = 1
                do j = decomp%zst(2),decomp%zen(2)
                    ii = 1
                    do i = decomp%zst(1),decomp%zen(1)
                        uexact(ii,jj,k) =  cos(x(i))*sin(z(k))*Ft
                        wexact(ii,jj,k) = -sin(x(i))*cos(z(k))*Ft
                        ii = ii + 1
                    end do 
                    jj = jj + 1
                end do 
            end do
            vexact = 0.d0 
        case(3)
            do k = decomp%zst(3),decomp%zen(3)
                jj = 1
                do j = decomp%zst(2),decomp%zen(2)
                    ii = 1
                    do i = decomp%zst(1),decomp%zen(1)
                        vexact(ii,jj,k) =  cos(y(j))*sin(z(k))*Ft  
                        wexact(ii,jj,k) = -sin(y(j))*cos(z(k))*Ft  
                        ii = ii + 1
                    end do 
                    jj = jj + 1
                end do 
            end do
            uexact = 0.d0 
        end select 

    end subroutine 

end module 

subroutine initfields_bolt(decomp, inputfile, delta_x, rho, ux, uy, uz)
    use kind_parameters, only:  rkind, clen
    use decomp_2d, only: decomp_info
    use constants, only: one, zero, two
    use d3q19_taylorgreen2D
    use exits, only: message
    use gridtools, only: linspace

    implicit none

    type(decomp_info), intent(in) :: decomp
    character(len=clen), intent(in) :: inputfile
    real(rkind), intent(in) :: delta_x
    real(rkind), dimension(:,:,:), intent(out) :: rho, ux, uy, uz  
    integer :: idx, i, j, k, ii, jj 

    nx = decomp%xsz(1)
    ny = decomp%ysz(2)
    nz = decomp%zsz(3)

    allocate(x(nx), y(ny), z(nz))

    x = linspace(0.d0,Lx-delta_x,nx)
    y = linspace(0.d0,Ly-delta_x,ny)
    z = linspace(0.d0,Ly-delta_x,nz)
    
    rho = one

    select case(TG_mode)
    case(1)
        do k = decomp%zst(3),decomp%zen(3)
            jj = 1
            do j = decomp%zst(2),decomp%zen(2)
                ii = 1
                do i = decomp%zst(1),decomp%zen(1)
                    ux(ii,jj,k) =  cos(x(i))*sin(y(j)) 
                    uy(ii,jj,k) = -sin(x(i))*cos(y(j)) 
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
                    ux(ii,jj,k) =  cos(x(i))*sin(z(k)) 
                    uz(ii,jj,k) = -sin(x(i))*cos(z(k)) 
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
                    uy(ii,jj,k) =  cos(y(j))*sin(z(k)) 
                    uz(ii,jj,k) = -sin(y(j))*cos(z(k)) 
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



