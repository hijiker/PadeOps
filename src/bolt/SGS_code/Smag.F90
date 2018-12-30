subroutine compute_tau_smag(this)
    class(d3q19), intent(inout) :: this
    integer :: i, j, k
    real(rkind) :: tmp, PiProd, SqrtTwoPiProd
    real(rkind) :: a, b, tau_hat_csq, xi4
    
    integer :: kst, ken

    if (this%isZPeriodic) then
        kst = 1
        ken = this%gp%zsz(3)
    else
        kst = 2
        ken = this%gp%zsz(3) - 1
    end if 

    select case (this%gradient_type)
    case (1)
        call this%compute_Pi()
        do k = kst,ken
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
                    
                    this%nusgs(i,j,k) = this%c_sgs_sq*sqrt(two*PiProd)*tmp
                    
                    this%tau(i,j,k) = (this%nusgs(i,j,k)+this%nu)*oneByCsq + half 
                end do 
            end do 
        end do
    case (2)
        call this%compute_duidxj(.false.)
        this%rbuffz1 = this%duidxj(:,:,:,1,1)*this%duidxj(:,:,:,1,1) + this%duidxj(:,:,:,2,2)*this%duidxj(:,:,:,2,2) &
                     + this%duidxj(:,:,:,3,3)*this%duidxj(:,:,:,3,3)
        
        this%rbuffz2 = half*(this%duidxj(:,:,:,1,2) + this%duidxj(:,:,:,2,1))
        this%rbuffz1 = this%rbuffz1 + two*this%rbuffz2*this%rbuffz2
        this%rbuffz2 = half*(this%duidxj(:,:,:,1,3) + this%duidxj(:,:,:,3,1))
        this%rbuffz1 = this%rbuffz1 + two*this%rbuffz2*this%rbuffz2
        this%rbuffz2 = half*(this%duidxj(:,:,:,2,3) + this%duidxj(:,:,:,3,2))
        this%rbuffz1 = this%rbuffz1 + two*this%rbuffz2*this%rbuffz2
        this%nusgs = this%c_sgs_sq*sqrt(two*this%rbuffz1)
        
        this%tau(:,:,kst:ken) = (this%nusgs(:,:,kst:ken) + this%nu)*oneByCsq + half
    
    case (3)
        call this%compute_Pi()
        tau_hat_csq = (this%nu/csq + half)*csq
        
        do k = kst,ken
            do j = 1,this%gp%zsz(2)
                ! !$omp simd
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
                    xi4 = -half*(-a + sqrt(a*a + b*this%c_sgs_sq))/(b*this%c_sgs_sq + 1.d-12)
                    this%nusgs(i,j,k) = xi4*xi4*this%c_sgs_sq*SqrtTwoPiProd
                    this%tau(i,j,k) = (this%nusgs(i,j,k)+this%nu)*oneByCsq + half 
                end do 
            end do 
        end do

    end select 

    if (.not. this%iszperiodic) then
        call getWall_nut(this%gp, this%delta_nu, this%ux, this%uy, this%uz, this%Re, this%tau_B, this%tau_T)
        this%tau(:,:,1) = this%tau_B
        this%tau(:,:,this%gp%zsz(3)) = this%tau_T
    end if 
end subroutine 
