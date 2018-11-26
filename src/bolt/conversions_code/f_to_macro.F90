    subroutine compute_macroscopic(this)
        use constants, only: half
        class(d3q19), intent(inout) :: this
        integer :: i, j, k
        real(rkind) :: onebyrho

        do k = 1,this%gp%zsz(3)
            do j = 1,this%gp%zsz(2)
                !$omp simd
                do i = 1,this%gp%zsz(1)
                    this%rho(i,j,k) = this%f(i,j,k,1 ) + this%f(i,j,k,2 ) + this%f(i,j,k,3 ) + this%f(i,j,k,4 ) &
                                 &  + this%f(i,j,k,5 ) + this%f(i,j,k,6 ) + this%f(i,j,k,7 ) + this%f(i,j,k,8 ) &
                                 &  + this%f(i,j,k,9 ) + this%f(i,j,k,10) + this%f(i,j,k,11) + this%f(i,j,k,12) &
                                 &  + this%f(i,j,k,13) + this%f(i,j,k,14) + this%f(i,j,k,15) + this%f(i,j,k,16) &
                                 &  + this%f(i,j,k,17) + this%f(i,j,k,18) + this%f(i,j,k,19)

                    onebyrho = one/(this%rho(i,j,k))

                    this%ux(i,j,k)  = (this%f(i,j,k,1 ) - this%f(i,j,k,2 ) + this%f(i,j,k,7 ) - this%f(i,j,k,8 ) &
                                 &  +  this%f(i,j,k,9 ) - this%f(i,j,k,10) + this%f(i,j,k,13) - this%f(i,j,k,14) &
                                 &  +  this%f(i,j,k,15) - this%f(i,j,k,16) + half*this%Fx)*onebyrho                           
                                                                                                               
                    this%uy(i,j,k)  = (this%f(i,j,k,3 ) - this%f(i,j,k,4 ) + this%f(i,j,k,7 ) - this%f(i,j,k,8 ) &
                                 &  +  this%f(i,j,k,11) - this%f(i,j,k,12) - this%f(i,j,k,13) + this%f(i,j,k,14) &
                                 &  +  this%f(i,j,k,17) - this%f(i,j,k,18) + half*this%Fy)*onebyrho                           
                                                                                                               
                    this%uz(i,j,k)  = (this%f(i,j,k,5 ) - this%f(i,j,k,6 ) + this%f(i,j,k,9 ) - this%f(i,j,k,10) &
                                 &  +  this%f(i,j,k,11) - this%f(i,j,k,12) - this%f(i,j,k,15) + this%f(i,j,k,16) &
                                 &  -  this%f(i,j,k,17) + this%f(i,j,k,18) + half*this%Fz)*onebyrho 
               
                end do 
            end do 
        end do 

    end subroutine 
