subroutine compute_tau_sigma(this)
    use constants, only: nine, six, three, pi
    class(d3q19), intent(inout) :: this
    real(rkind) :: G11, G12, G13, G22, G23, G33
    real(rkind) :: I1, I2, I3, I1sq, I1cu
    real(rkind) :: alpha1, alpha2, alpha1sqrt
    real(rkind) :: alpha1tmp, alpha3
    real(rkind) :: sigma1, sigma2, sigma3, sigma1sq, scale_fact
    integer :: i,j,k, kst, ken

    call this%compute_duidxj(.true.)
    
    scale_fact = ((this%c_sgs*this%delta_x)**2)/this%delta_nu
    if (this%isZPeriodic) then
        kst = 1
        ken = this%gp%zsz(3)
    else
        kst = 2
        ken = this%gp%zsz(3) - 1
    end if
    
    do k = kst,ken
       do j = 1,this%gp%zsz(2)
          !$omp simd 
          do i = 1,this%gp%zsz(1)
             
             G11 = this%duidxj(i,j,k,1,1)*this%duidxj(i,j,k,1,1) + this%duidxj(i,j,k,2,1)*this%duidxj(i,j,k,2,1) + this%duidxj(i,j,k,3,1)*this%duidxj(i,j,k,3,1)
             G12 = this%duidxj(i,j,k,1,1)*this%duidxj(i,j,k,1,2) + this%duidxj(i,j,k,2,1)*this%duidxj(i,j,k,2,2) + this%duidxj(i,j,k,3,1)*this%duidxj(i,j,k,3,2)
             G13 = this%duidxj(i,j,k,1,1)*this%duidxj(i,j,k,1,3) + this%duidxj(i,j,k,2,1)*this%duidxj(i,j,k,2,3) + this%duidxj(i,j,k,3,1)*this%duidxj(i,j,k,3,3)
             G22 = this%duidxj(i,j,k,1,2)*this%duidxj(i,j,k,1,2) + this%duidxj(i,j,k,2,2)*this%duidxj(i,j,k,2,2) + this%duidxj(i,j,k,3,2)*this%duidxj(i,j,k,3,2)
             G23 = this%duidxj(i,j,k,1,2)*this%duidxj(i,j,k,1,3) + this%duidxj(i,j,k,2,2)*this%duidxj(i,j,k,2,3) + this%duidxj(i,j,k,3,2)*this%duidxj(i,j,k,3,3)
             G33 = this%duidxj(i,j,k,1,3)*this%duidxj(i,j,k,1,3) + this%duidxj(i,j,k,2,3)*this%duidxj(i,j,k,2,3) + this%duidxj(i,j,k,3,3)*this%duidxj(i,j,k,3,3)


             I1   = G11 + G22 + G33
             I1sq = I1*I1
             I1cu = I1sq*I1

             I2 = -G11*G11 - G22*G22 - G33*G33
             I2 = I2 - two*G12*G12 - two*G13*G13
             I2 = I2 - two*G23*G23
             I2 = I2 + I1sq
             I2 = half*I2

             I3 = G11*(G22*G33 - G23*G23)
             I3 = I3 + G12*(G13*G23 - G12*G33)
             I3 = I3 + G13*(G12*G23 - G22*G13)

             alpha1 = I1sq/nine - I2/three
             alpha1 = max(alpha1,zero)

             alpha2 = I1cu/27._rkind - I1*I2/six + I3/two
             alpha1sqrt = sqrt(alpha1)
             alpha1tmp = alpha1*alpha1sqrt
             alpha1tmp = alpha2/(alpha1tmp + 1.d-13)
             alpha1tmp = min(alpha1tmp,one)
             alpha1tmp = max(alpha1tmp,-one)
             alpha1tmp = acos(alpha1tmp)
             alpha3 = (one/three)*(alpha1tmp)

             sigma1sq = I1/three + two*alpha1sqrt*cos(alpha3)
             sigma1sq = max(sigma1sq,zero)
             sigma1 = sqrt(sigma1sq)

             sigma2 = pi/three + alpha3
             sigma2 = (-two)*alpha1sqrt*cos(sigma2)
             sigma2 = sigma2 + I1/three
             sigma2 = max(sigma2,zero)
             sigma2 = sqrt(sigma2)
             
             sigma3 = pi/three - alpha3
             sigma3 = (-two)*alpha1sqrt*cos(sigma3)
             sigma3 = sigma3 + I1/three
             sigma3 = max(sigma3,zero)
             sigma3 = sqrt(sigma3)

             this%nusgs(i,j,k) = (sigma3*(sigma1 - sigma2)*(sigma2 - sigma3)/(sigma1sq + 1.d-15))*scale_fact
             this%tau(i,j,k) = (this%nusgs(i,j,k)+this%nu)*oneByCsq + half 

          end do 
       end do 
    end do

    if (.not. this%iszperiodic) then
        call getWall_nut(this%gp, this%delta_nu, this%ux, this%uy, this%uz, this%Re, this%tau_B, this%tau_T)
        this%tau(:,:,1) = this%tau_B
        this%tau(:,:,this%gp%zsz(3)) = this%tau_T
    end if 



end subroutine 
