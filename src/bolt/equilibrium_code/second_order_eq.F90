    pure subroutine get_Feq_2ndOrder(ux,uy,uz,rho,fidx,Qtensor,feq)
        real(rkind), intent(in) :: ux, uy, uz, rho
        integer, intent(in) :: fidx
        real(rkind), dimension(3,3), intent(in) :: Qtensor
        real(rkind), intent(out) :: feq

        real(rkind) :: first, second

        first = cx(fidx)*ux + cy(fidx)*uy + cz(fidx)*uz

        second =     ux*ux*Qtensor(1,1) +     uy*uy*Qtensor(2,2) +     uz*uz*Qtensor(3,3) &
             & + two*ux*uy*Qtensor(1,2) + two*ux*uz*Qtensor(1,3) + two*uy*uz*Qtensor(2,3) 
        
        feq = w(fidx)*rho*(one + onebycsq*first + oneby2c4*second)

    end subroutine 

    pure subroutine get_ForceSource_2ndOrder(ux,uy,uz,ForceX,ForceY,ForceZ,fidx,Qtensor,Fsource)
        real(rkind), intent(in) :: ux, uy, uz, ForceX, ForceY, ForceZ
        integer, intent(in) :: fidx
        real(rkind), dimension(3,3), intent(in) :: Qtensor
        real(rkind), intent(out) :: Fsource

        real(rkind) :: first, second

        first = cx(fidx)*ForceX + cy(fidx)*ForceY + cz(fidx)*ForceZ

        second = (Qtensor(1,1)*ux + Qtensor(1,2)*uy + Qtensor(1,3)*uz)*ForceX &
               + (Qtensor(2,1)*ux + Qtensor(2,2)*uy + Qtensor(2,3)*uz)*ForceY &
               + (Qtensor(3,1)*ux + Qtensor(3,2)*uy + Qtensor(3,3)*uz)*ForceZ
      
        Fsource = w(fidx)*(onebycsq*first + oneby2c4*second)

    end subroutine 
    
    pure subroutine get_ForceSource_Guo(ux,uy,uz,ForceX,ForceY,ForceZ,fidx,Fsource)
        real(rkind), intent(in) :: ux, uy, uz, ForceX, ForceY, ForceZ
        integer, intent(in) :: fidx
        real(rkind), intent(out) :: Fsource

        real(rkind) :: first, second, cdotu

        first = (cx(fidx) - ux)*ForceX + (cy(fidx) - uy)*ForceY + (cz(fidx) - uz)*ForceZ
        cdotu = (cx(fidx)*ux + cy(fidx)*uy + cz(fidx)*uz) 
        second = cdotu*(cx(fidx)*ForceX + cy(fidx)*ForceY + cz(fidx)*ForceZ)

        Fsource = w(fidx)*(onebycsq*first + onebycsq*onebycsq*second)

    end subroutine 
