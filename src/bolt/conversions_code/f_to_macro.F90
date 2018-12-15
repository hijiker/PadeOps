    subroutine compute_macroscopic(this)
        use constants, only: half
        class(d3q19), intent(inout) :: this
        integer :: k, kst, ken

        if (this%isZperiodic) then
            kst = 1
            ken = this%gp%zsz(3)
        else
            kst = 2
            ken = this%gp%zsz(3) - 1
        end if 

        do k = 1,this%gp%zsz(3)
            call this%compute_macroscopic_hplane(k)
        end do 

    end subroutine 
   
    subroutine compute_macroscopic_hplane(this, k)
        class(d3q19), intent(inout) :: this
        integer, intent(in) :: k
        integer :: i, j
        real(rkind) :: onebyrho
            
        do j = 1,this%gp%zsz(2)
            !$omp simd
            do i = 1,this%gp%zsz(1)
                this%rho(i,j,k) = this%f(i,j,k,velOrder(1 )) + this%f(i,j,k,velOrder(2 )) + this%f(i,j,k,velOrder(3 )) + this%f(i,j,k,velOrder(4 )) &
                             &  + this%f(i,j,k,velOrder(5 )) + this%f(i,j,k,velOrder(6 )) + this%f(i,j,k,velOrder(7 )) + this%f(i,j,k,velOrder(8 )) &
                             &  + this%f(i,j,k,velOrder(9 )) + this%f(i,j,k,velOrder(10)) + this%f(i,j,k,velOrder(11)) + this%f(i,j,k,velOrder(12)) &
                             &  + this%f(i,j,k,velOrder(13)) + this%f(i,j,k,velOrder(14)) + this%f(i,j,k,velOrder(15)) + this%f(i,j,k,velOrder(16)) &
                             &  + this%f(i,j,k,velOrder(17)) + this%f(i,j,k,velOrder(18)) + this%f(i,j,k,velOrder(19))

                onebyrho = one/(this%rho(i,j,k))

                this%ux(i,j,k)  = (this%f(i,j,k,velOrder(1 )) - this%f(i,j,k,velOrder(2 )) + this%f(i,j,k,velOrder(7 )) - this%f(i,j,k,velOrder(8 )) &
                             &  +  this%f(i,j,k,velOrder(9 )) - this%f(i,j,k,velOrder(10)) + this%f(i,j,k,velOrder(13)) - this%f(i,j,k,velOrder(14)) &
                             &  +  this%f(i,j,k,velOrder(15)) - this%f(i,j,k,velOrder(16)) + half*this%Fx(i,j,k))*onebyrho                           
                                                                                                           
                this%uy(i,j,k)  = (this%f(i,j,k,velOrder(3 )) - this%f(i,j,k,velOrder(4 )) + this%f(i,j,k,velOrder(7 )) - this%f(i,j,k,velOrder(8 )) &
                             &  +  this%f(i,j,k,velOrder(11)) - this%f(i,j,k,velOrder(12)) - this%f(i,j,k,velOrder(13)) + this%f(i,j,k,velOrder(14)) &
                             &  +  this%f(i,j,k,velOrder(17)) - this%f(i,j,k,velOrder(18)) + half*this%Fy(i,j,k))*onebyrho                           
                                                                                                           
                this%uz(i,j,k)  = (this%f(i,j,k,velOrder(5 )) - this%f(i,j,k,velOrder(6 )) + this%f(i,j,k,velOrder(9 )) - this%f(i,j,k,velorder(10)) &
                             &  +  this%f(i,j,k,velOrder(11)) - this%f(i,j,k,velOrder(12)) - this%f(i,j,k,velOrder(15)) + this%f(i,j,k,velorder(16)) &
                             &  -  this%f(i,j,k,velOrder(17)) + this%f(i,j,k,velOrder(18)) + half*this%Fz(i,j,k))*onebyrho 
           
            end do 
        end do 

    end subroutine 
