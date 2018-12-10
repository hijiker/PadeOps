subroutine collision_MRT_Force(this)
    use MRT_D3Q19_coefficients, only: get_MRT_lambda 
    class(d3q19), intent(inout) :: this
    real(rkind), dimension(nvels) :: mvec, ftmp, lambda_vec, meq, Force, Mf
    integer :: i, j, k, idx

    do k = 1,this%gp%zsz(3)
        do j = 1,this%gp%zsz(2)
            do i = 1,this%gp%zsz(1)
                ftmp(2:nvels) = this%f(i,j,k,1:nvels-1)
                ftmp(1) = this%f(i,j,k,nvels)
                !ftmp = this%f(i,j,k,:)

                call fspace_to_mspace(ftmp,mvec)
                
                call get_ForceSource_Guo(this%ux(i,j,k),this%uy(i,j,k),this%uz(i,j,k), &
                        & this%Fx(i,j,k), this%Fy(i,j,k), this%Fz(i,j,k), nvels,  Force(1))
                !$omp simd
                do idx = 1,nvels-1 
                        call get_ForceSource_Guo(this%ux(i,j,k),this%uy(i,j,k),this%uz(i,j,k), &
                                & this%Fx(i,j,k), this%Fy(i,j,k), this%Fz(i,j,k), idx, Force(idx+1))
                end do 
                call transform_Force_to_mspace(Force,Mf)
                
                call get_MRT_lambda(this%tau(i,j,k), lambda_vec)
                call compute_Meq(this%rho(i,j,k),this%ux(i,j,k),this%uy(i,j,k),this%uz(i,j,k),meq)

                mvec = (one - lambda_vec)*mvec + lambda_vec*meq + (one - half*lambda_vec)*Mf
                
                call mspace_to_fspace(mvec,ftmp)
                
                !this%f(i,j,k,:) = ftmp
                this%f(i,j,k,1:nvels-1) = ftmp(2:nvels)
                this%f(i,j,k,nvels) = ftmp(1)
            end do
        end do 
    end do


end subroutine

pure subroutine fspace_to_mspace(f,mvec)
    use MRT_D3Q19_coefficients, only: M
    real(rkind), dimension(nvels), intent(in) :: f
    real(rkind), dimension(nvels), intent(out) :: mvec
    integer :: i
    
    mvec = zero 
    do i = 1,nvels
        mvec = mvec+ f(i)*M(:,i) 
    end do 
end subroutine 

pure subroutine mspace_to_fspace(mvec,f)
    use MRT_D3Q19_coefficients, only: Minv
    real(rkind), dimension(nvels), intent(in) :: mvec 
    real(rkind), dimension(nvels), intent(out) :: f
    integer :: i

    f = zero
    do i = 1,nvels
        f = f + mvec(i)*Minv(:,i) 
    end do

end subroutine 

pure subroutine transform_Force_to_mspace(Force,Mf)
    use MRT_D3Q19_coefficients, only: M
    real(rkind), dimension(nvels), intent(in) :: Force
    real(rkind), dimension(nvels), intent(out) :: Mf
    integer :: i
    
    Mf = zero 
    do i = 1,nvels
        Mf = Mf + Force(i)*M(:,i) 
    end do 

end subroutine 


pure subroutine compute_Meq(rho,ux,uy,uz,meq)
    real(rkind), intent(in) :: rho, ux, uy, uz
    real(rkind), dimension(nvels), intent(out) :: meq
    
    meq(1 ) = rho
    meq(2 ) = -11*rho + 19*rho*(ux*ux + uy*uy + uz*uz)
    meq(3 ) = 3*rho - (11.d0/2.d0)*rho*(ux*ux + uy*uy + uz*uz)
    meq(4 ) = rho*ux
    meq(5 ) = -(2.d0/3.d0)*rho*ux
    meq(6 ) = rho*uy
    meq(7 ) = -(2.d0/3.d0)*rho*uy
    meq(8 ) = rho*uz
    meq(9 ) = -(2.d0/3.d0)*rho*uz
    meq(10) = 2*rho*ux*ux - rho*uy*uy - rho*uz*uz
    meq(11) = -rho*ux*ux + 0.5d0*uy*uy + 0.5d0*uz*uz
    meq(12) = rho*uy*uy - rho*uz*uz
    meq(13) = -0.5d0*rho*uy*uy + 0.5d0*rho*uz*uz
    meq(14) = rho*ux*uy
    meq(15) = rho*uz*uy
    meq(16) = rho*ux*uz
    meq(17) = zero
    meq(18) = zero
    meq(19) = zero

end subroutine
