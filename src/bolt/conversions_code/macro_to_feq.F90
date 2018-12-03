    subroutine initialize_f_to_feq(this)
        class(d3q19), intent(inout) :: this
        integer :: i, j, k, idx 
        real(rkind) :: uin, vin, win

        do idx = 1,nvels
            do k = 1,this%gp%zsz(3)
                do j = 1,this%gp%zsz(2)
                    !$omp simd 
                    do i = 1,this%gp%zsz(1)
                        uin = this%ux(i,j,k) - this%Fx(i,j,k)/(two*this%rho(i,j,k))
                        vin = this%uy(i,j,k) - this%Fy(i,j,k)/(two*this%rho(i,j,k))
                        win = this%uz(i,j,k) - this%Fz(i,j,k)/(two*this%rho(i,j,k))
                        call get_Feq_2ndOrder(uin,vin,win, this%rho(i,j,k), &
                            & idx,this%Qtensor(:,:,idx),this%f(i,j,k,idx))
                        
                    end do 
                end do 
            end do 
        end do 

    end subroutine
