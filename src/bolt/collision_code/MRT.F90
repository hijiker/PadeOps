subroutine collision_MRT_Force(this)
    use MRT_D3Q19_coefficients, only: get_MRT_lambda 
    class(d3q19), intent(inout) :: this
    integer :: i, j, k, idx

    do k = 1,this%gp%zsz(3)
        do j = 1,this%gp%zsz(2)
            do i = 1,this%gp%zsz(1)
                this%f_mrt = this%f(i,j,k,:)

                call this%fspace_to_mspace()
                
                !$omp simd
                do idx = 1,nvels 
                        call get_ForceSource_Guo(this%ux(i,j,k),this%uy(i,j,k),this%uz(i,j,k), &
                                & this%Fx(i,j,k), this%Fy(i,j,k), this%Fz(i,j,k), idx, this%force_mrt(idx))
                end do
 
                call this%transform_Force_to_mspace()
                
                call get_MRT_lambda(this%tau(i,j,k), this%lambda_mrt)
                
                call this%compute_Meq(this%rho(i,j,k),this%ux(i,j,k),this%uy(i,j,k),this%uz(i,j,k))

                this%m_mrt = (one - this%lambda_mrt)*this%m_mrt + this%lambda_mrt*this%meq_mrt + (one - half*this%lambda_mrt)*this%Mf_mrt
                
                call this%mspace_to_fspace()
                
                this%f(i,j,k,:) = this%f_mrt
            end do
        end do 
    end do


end subroutine

pure subroutine fspace_to_mspace(this)
    use MRT_D3Q19_coefficients, only: Mmat
    class(d3q19), intent(inout) :: this
    integer :: i
    
    this%m_mrt = zero 
    do i = 1,nvels
        this%m_mrt = this%m_mrt+ this%f_mrt(i)*Mmat(:,i) 
    end do 
end subroutine 

pure subroutine mspace_to_fspace(this)
    use MRT_D3Q19_coefficients, only: Minv
    class(d3q19), intent(inout) :: this
    integer :: i

    this%f_mrt = zero
    do i = 1,nvels
        this%f_mrt = this%f_mrt + this%m_mrt(i)*Minv(:,i) 
    end do

end subroutine 

pure subroutine transform_Force_to_mspace(this)
    use MRT_D3Q19_coefficients, only: Mmat
    class(d3q19), intent(inout) :: this
    integer :: i
    
    this%Mf_mrt = zero 
    do i = 1,nvels
        this%Mf_mrt = this%Mf_mrt + this%force_mrt(i)*Mmat(:,i) 
    end do 

end subroutine 


pure subroutine compute_Meq(this,rho,ux,uy,uz)
    use MRT_D3Q19_coefficients, only: Mmat
    class(d3q19), intent(inout) :: this
    real(rkind), intent(in) :: rho, ux, uy, uz
    !real(rkind), dimension(nvels) :: feq
    !integer :: idx, i

    !!$omp simd
    !do idx = 1,nvels
    !    call get_Feq_2ndOrder(ux,uy,uz,rho,idx,this%Qtensor(:,:,idx),feq(idx))
    !end do 
    !
    !this%meq_mrt = zero  
    !do i = 1,nvels
    !    this%meq_mrt = this%meq_mrt + feq(i)*Mmat(:,i) 
    !end do 

    this%meq_mrt(1 ) = rho
    this%meq_mrt(2 ) = -11*rho + 19*rho*(ux*ux + uy*uy + uz*uz)
    this%meq_mrt(3 ) = 3*rho - (11.d0/2.d0)*rho*(ux*ux + uy*uy + uz*uz)
    this%meq_mrt(4 ) = rho*ux
    this%meq_mrt(5 ) = -(2.d0/3.d0)*rho*ux
    this%meq_mrt(6 ) = rho*uy
    this%meq_mrt(7 ) = -(2.d0/3.d0)*rho*uy
    this%meq_mrt(8 ) = rho*uz
    this%meq_mrt(9 ) = -(2.d0/3.d0)*rho*uz
    this%meq_mrt(10) = 2*rho*ux*ux - rho*uy*uy - rho*uz*uz
    this%meq_mrt(11) = -rho*ux*ux + 0.5d0*uy*uy + 0.5d0*uz*uz
    this%meq_mrt(12) = rho*uy*uy - rho*uz*uz
    this%meq_mrt(13) = -0.5d0*rho*uy*uy + 0.5d0*rho*uz*uz
    this%meq_mrt(14) = rho*ux*uy
    this%meq_mrt(15) = rho*uz*uy
    this%meq_mrt(16) = rho*ux*uz
    this%meq_mrt(17) = zero
    this%meq_mrt(18) = zero
    this%meq_mrt(19) = zero

end subroutine
