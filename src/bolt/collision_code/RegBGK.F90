    subroutine RegBGK_Force(this)
        class(d3q19), intent(inout) :: this
        real(rkind), dimension(nvels) :: feq, Fvals
        integer :: i, j, k, idx
        real(rkind), dimension(6) :: Preg
        real(rkind) :: fneq, freg, OneByTau

        do k = 1,this%gp%zsz(3)
            do j = 1,this%gp%zsz(2)
                do i = 1,this%gp%zsz(1)
                    Preg = zero
                    do idx = 1,nvels
                        call get_Feq_2ndOrder(this%ux(i,j,k),this%uy(i,j,k),this%uz(i,j,k), &
                                & this%rho(i,j,k),idx,this%Qtensor(:,:,idx),feq(idx))
                        
                        fneq = this%f(i,j,k,idx) - feq(idx)
                        
                        Preg(1) = Preg(1) + cx(idx)*cx(idx)*(fneq)
                        Preg(2) = Preg(2) + cx(idx)*cy(idx)*(fneq)
                        Preg(3) = Preg(3) + cx(idx)*cz(idx)*(fneq)
                        Preg(4) = Preg(4) + cy(idx)*cy(idx)*(fneq)
                        Preg(5) = Preg(5) + cy(idx)*cz(idx)*(fneq)
                        Preg(6) = Preg(6) + cz(idx)*cz(idx)*(fneq)
                        
                    end do
                   
                    oneByTau = one/this%tau(i,j,k)
                    do idx = 1,nvels
                        freg = w(idx)*oneby2c4*(Preg(1)*this%Qtensor(1,1,idx) + two*Preg(2)*this%Qtensor(1,2,idx) + &
                                  & two*Preg(3)*this%Qtensor(1,3,idx) + Preg(4)*this%Qtensor(2,2,idx)  +            &
                                  & two*Preg(5)*this%Qtensor(2,3,idx) + Preg(6)*this%Qtensor(3,3,idx)) -            &
                                  & w(idx)*half*onebycsq*(this%Fx*cx(idx) + this%Fy*cy(idx) + this%Fz*cz(idx)) 
                   
                        ! Collision
                        this%f(i,j,k,idx) = feq(idx) + (one - oneByTau)*freg + (one - half*oneBytau)*Fvals(idx)
                    end do 
                end do
            end do 
        end do 
   
    end subroutine 

