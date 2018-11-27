    subroutine collision_BGK_Force(this)
        class(d3q19), intent(inout) :: this
        integer :: i, j, k, idx 
        real(rkind) :: OneByTau, feq, Force

        do idx = 1,nvels
            do k = 1,this%gp%zsz(3)
                do j = 1,this%gp%zsz(2)
                    !$omp simd 
                    do i = 1,this%gp%zsz(1)
                        call get_Feq_2ndOrder(this%ux(i,j,k),this%uy(i,j,k),this%uz(i,j,k), &
                                & this%rho(i,j,k),idx,this%Qtensor(:,:,idx),feq)
                        
                        call get_ForceSource_2ndOrder(this%ux(i,j,k),this%uy(i,j,k),this%uz(i,j,k), &
                                & this%Fx, this%Fy, this%Fz, idx, &
                                & this%Qtensor(:,:,idx), Force)
                        
                        oneByTau = one/this%tau(i,j,k)
                       
                        this%f(i,j,k,idx) = (one - oneByTau)*this%f(i,j,k,idx) + oneByTau*feq &
                                        & + (one - half*oneBytau)*Force
                    end do 
                end do 
            end do 
        end do 
        
    end subroutine


    subroutine collision_BGK(this)
        class(d3q19), intent(inout) :: this
        integer :: i, j, k, idx 
        real(rkind) :: OneByTau, feq, Force

        do idx = 1,nvels
            do k = 1,this%gp%zsz(3)
                do j = 1,this%gp%zsz(2)
                    !$omp simd 
                    do i = 1,this%gp%zsz(1)
                        call get_Feq_2ndOrder(this%ux(i,j,k),this%uy(i,j,k),this%uz(i,j,k), &
                                & this%rho(i,j,k),idx,this%Qtensor(:,:,idx),feq)
                        
                        call get_ForceSource_2ndOrder(this%ux(i,j,k),this%uy(i,j,k),this%uz(i,j,k), &
                                & this%Fx, this%Fy, this%Fz, idx, &
                                & this%Qtensor(:,:,idx), Force)
                        
                        oneByTau = one/this%tau(i,j,k)
                        this%f(i,j,k,idx) = (one - oneByTau)*this%f(i,j,k,idx) + oneByTau*feq
                    end do 
                end do 
            end do 
        end do 

    end subroutine
