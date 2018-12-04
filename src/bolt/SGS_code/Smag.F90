subroutine compute_tau_smag(this)
    class(d3q19), intent(inout) :: this
    integer :: i, j, k
    real(rkind) :: tmp, PiProd

    call this%compute_Pi()
    do k = 1,this%gp%zsz(3)
        do j = 1,this%gp%zsz(2)
            !$omp simd
            do i = 1,this%gp%zsz(1)
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

    call getWall_nut(this%gp, this%delta_nu, this%ux, this%uy, this%uz, this%Re, this%tau_B, this%tau_T)
    this%tau(:,:,1) = this%tau_B
    this%tau(:,:,this%gp%zsz(3)) = this%tau_T

end subroutine 
