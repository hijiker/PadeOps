subroutine compute_tau_smag(this)
    use reductions 
    class(d3q19), intent(inout) :: this
    integer :: i, j, k
    real(rkind) :: tmp, PiProd, SqrtTwoPiProd
    real(rkind) :: a, b, tau_hat_csq, xi4
    select case (this%gradient_type)
    case (1)
        call this%compute_Pi()
        do k = 1,this%gp%zsz(3)
            do j = 1,this%gp%zsz(2)
                !$omp simd
                do i = 1,this%gp%zsz(1)
                    ! use previous nu_sgs
                    tmp = one/(this%tau(i,j,k)*csq)
                    PiProd = this%PiTensor(1,i,j,k)*this%PiTensor(1,i,j,k) &
                           + this%PiTensor(4,i,j,k)*this%PiTensor(4,i,j,k) &
                           + this%PiTensor(6,i,j,k)*this%PiTensor(6,i,j,k) &
                           + two*this%PiTensor(2,i,j,k)*this%PiTensor(2,i,j,k) &
                           + two*this%PiTensor(3,i,j,k)*this%PiTensor(3,i,j,k) &
                           + two*this%PiTensor(5,i,j,k)*this%PiTensor(5,i,j,k)
                    
                    this%nusgs(i,j,k) = C_Smag_sq*sqrt(two*PiProd)*tmp
                    
                    this%tau(i,j,k) = (this%nusgs(i,j,k)+this%nu)*oneByCsq + half 
                end do 
            end do 
        end do
    case (2)
        call this%compute_duidxj()
        this%rbuffz1 = this%duidxj(:,:,:,1,1)*this%duidxj(:,:,:,1,1) + this%duidxj(:,:,:,2,2)*this%duidxj(:,:,:,2,2) &
                     + this%duidxj(:,:,:,3,3)*this%duidxj(:,:,:,3,3)
        
        this%rbuffz2 = half*(this%duidxj(:,:,:,1,2) + this%duidxj(:,:,:,2,1))
        this%rbuffz1 = this%rbuffz1 + two*this%rbuffz2*this%rbuffz2
        this%rbuffz2 = half*(this%duidxj(:,:,:,1,3) + this%duidxj(:,:,:,3,1))
        this%rbuffz1 = this%rbuffz1 + two*this%rbuffz2*this%rbuffz2
        this%rbuffz2 = half*(this%duidxj(:,:,:,2,3) + this%duidxj(:,:,:,3,2))
        this%rbuffz1 = this%rbuffz1 + two*this%rbuffz2*this%rbuffz2
        this%nusgs = C_smag_sq*sqrt(two*this%rbuffz1)
        
        this%tau = (this%nusgs + this%nu)*oneByCsq + half
    
    case (3)
        call this%compute_Pi()
        tau_hat_csq = (this%nu/csq + half)*csq
        
        do k = 1,this%gp%zsz(3)
            do j = 1,this%gp%zsz(2)
                !$omp simd
                do i = 1,this%gp%zsz(1)
                    PiProd = this%PiTensor(1,i,j,k)*this%PiTensor(1,i,j,k) &
                           + this%PiTensor(4,i,j,k)*this%PiTensor(4,i,j,k) &
                           + this%PiTensor(6,i,j,k)*this%PiTensor(6,i,j,k) &
                           + two*this%PiTensor(2,i,j,k)*this%PiTensor(2,i,j,k) &
                           + two*this%PiTensor(3,i,j,k)*this%PiTensor(3,i,j,k) &
                           + two*this%PiTensor(5,i,j,k)*this%PiTensor(5,i,j,k)
                    
                    SqrtTwoPiProd = sqrt(2.d0*PiProd)

                    a = this%rho(i,j,k)*tau_hat_csq
                    b = 2.d0*SqrtTwoPiProd*this%rho(i,j,k)

                    xi4 = -half*(-a + sqrt(a*a + b*C_smag_sq))/(b*C_smag_sq)
                    this%nusgs(i,j,k) = xi4*xi4*C_smag_sq*SqrtTwoPiProd
                    this%tau(i,j,k) = (this%nusgs(i,j,k)+this%nu)*oneByCsq + half 
                end do 
            end do 
        end do

    end select 

    call getWall_nut(this%gp, this%delta_nu, this%ux, this%uy, this%uz, this%Re, this%tau_B, this%tau_T)
    this%tau(:,:,1) = this%tau_B
    this%tau(:,:,this%gp%zsz(3)) = this%tau_T

end subroutine 
    
subroutine compute_duidxj(this)
    class(d3q19), intent(inout) :: this

    call this%get_gradient(this%ux,this%duidxj(:,:,:,1,1),this%duidxj(:,:,:,1,2),this%duidxj(:,:,:,1,3))
    call this%get_gradient(this%uy,this%duidxj(:,:,:,2,1),this%duidxj(:,:,:,2,2),this%duidxj(:,:,:,2,3))
    call this%get_gradient(this%uz,this%duidxj(:,:,:,3,1),this%duidxj(:,:,:,3,2),this%duidxj(:,:,:,3,3))
   
end subroutine

subroutine get_gradient(this, u, dudx, dudy, dudz)
    class(d3q19), intent(inout) :: this
    real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)), intent(in ) :: u
    real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)), intent(out) :: dudx, dudy, dudz

    call this%get_ddz(u,dudz)

    call transpose_z_to_y(u,this%rbuffy1,this%gp)
    call this%get_ddy(this%rbuffy1,this%rbuffy2)
    call transpose_y_to_z(this%rbuffy2,dudy)

    call transpose_y_to_x(this%rbuffy1,this%rbuffx1)
    call this%get_ddx(this%rbuffx1,this%rbuffx2)
    call transpose_x_to_y(this%rbuffx2,this%rbuffy2,this%gp)
    call transpose_y_to_z(this%rbuffy2,dudx,this%gp)

end subroutine 


subroutine get_ddz(this,u,dudz)
    class(d3q19), intent(in) :: this
    real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)), intent(in ) :: u
    real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)), intent(out) :: dudz
    integer :: nz

    nz = this%gp%zsz(3)
    
    ! dudz
    dudz(:,:,2:nz-1) = half*(u(:,:,3:nz) - u(:,:,1:nz-2))
    dudz(:,:,1) = (-3.d0/2.d0)*u(:,:,1) + 2.d0*u(:,:,2) - 0.5d0*(u(:,:,3)) 
    dudz(:,:,nz) = -((-3.d0/2.d0)*u(:,:,nz) + 2.d0*u(:,:,nz-1) - 0.5d0*(u(:,:,nz-2)))
   
end subroutine


subroutine get_ddy(this,u,dudy)
    class(d3q19), intent(in) :: this
    real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3)), intent(in ) :: u
    real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3)), intent(out) :: dudy
    integer :: ny 

    ny = this%gp%ysz(2)

    dudy(:,2:ny-1,:) = half*(u(:,3:ny,:) - u(:,1:ny-2,:))
    dudy(:,1,:) = half*(u(:,2,:) - u(:,ny,:))
    dudy(:,ny,:) = half*(u(:,1,:) - u(:,ny-1,:))
    
end subroutine

subroutine get_ddx(this,u,dudx)
    class(d3q19), intent(in) :: this
    real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in ) :: u
    real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: dudx
    integer :: nx 

    nx = this%gp%xsz(1)

    dudx(2:nx-1,:,:) = half*(u(3:nx,:,:) - u(1:nx-2,:,:))
    dudx(1,:,:) = half*(u(2,:,:) - u(nx,:,:))
    dudx(nx,:,:) = half*(u(1,:,:) - u(nx-1,:,:))
end subroutine
