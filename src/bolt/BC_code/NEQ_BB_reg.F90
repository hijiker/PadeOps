    subroutine NEQ_BB_Reg(this) 
        class(d3q19), intent(inout) :: this
        integer :: vid

        ! Impose BOTTOM BC
        call this%compute_rhoBC(fminus,this%wBot,this%f(:,:,1,:),-1)
        call this%compute_f_BC(this%uBot, this%vBot, this%wBot,this%f(:,:,1,:))
        do vid = 1,size(fplus)
            this%fneqBC(:,:,fplus(vid)) = this%fneqBC(:,:,fminus(vid))
        end do 
        call this%RegularizeFneq_BC()
        this%f(:,:,1,:) = this%feqBC + this%fneqBC
        
        ! Impose TOP BC
        call this%compute_rhoBC(fplus,this%wBot,this%f(:,:,this%nz,:),1)
        call this%compute_f_BC(this%uBot, this%vBot, this%wBot,this%f(:,:,this%nz,:))
        do vid = 1,size(fminus)
            this%fneqBC(:,:,fminus(vid)) = this%fneqBC(:,:,fplus(vid))
        end do 
        call this%RegularizeFneq_BC()
        this%f(:,:,this%nz,:) = this%feqBC + this%fneqBC

    end subroutine 

    subroutine RegularizeFneq_BC(this)
        class(d3q19), intent(inout) :: this
        integer :: vid, i, j

        this%PiBC = zero 
        do vid = 1,nvels
            do j = 1,this%gp%zsz(2)
                !$omp simd
                do i = 1,this%gp%zsz(1)
                    this%PiBC(i,j,1,1) = this%PiBC(i,j,1,1) + cx(vid)*cx(vid)*this%fneqBC(i,j,vid)
                    this%PiBC(i,j,1,2) = this%PiBC(i,j,1,2) + cx(vid)*cy(vid)*this%fneqBC(i,j,vid)
                    this%PiBC(i,j,1,3) = this%PiBC(i,j,1,3) + cx(vid)*cz(vid)*this%fneqBC(i,j,vid)
                    
                    !this%PiBC(i,j,2,1) = this%PiBC(i,j,2,1) + cy(vid)*cx(vid)*this%fneqBC(i,j,vid)
                    this%PiBC(i,j,2,2) = this%PiBC(i,j,2,2) + cy(vid)*cy(vid)*this%fneqBC(i,j,vid)
                    this%PiBC(i,j,2,3) = this%PiBC(i,j,2,3) + cy(vid)*cz(vid)*this%fneqBC(i,j,vid)
                    
                    !this%PiBC(i,j,3,1) = this%PiBC(i,j,3,1) + cz(vid)*cx(vid)*this%fneqBC(i,j,vid)
                    !this%PiBC(i,j,3,2) = this%PiBC(i,j,3,2) + cz(vid)*cy(vid)*this%fneqBC(i,j,vid)
                    this%PiBC(i,j,3,3) = this%PiBC(i,j,3,3) + cz(vid)*cz(vid)*this%fneqBC(i,j,vid)
                end do 
            end do 
        end do 

        do vid = 1,nvels
            !this%fneqBC(:,:,vid) = (w(vid)*oneby2c4)*(this%PiBC(:,:,1,1)*this%Qtensor(1,1,vid) &
            !                     & + this%PiBC(:,:,1,2)*this%Qtensor(1,2,vid) + this%PiBC(:,:,1,3)*this%Qtensor(1,3,vid) & 
            !                     & + this%PiBC(:,:,2,1)*this%Qtensor(2,1,vid) + this%PiBC(:,:,2,2)*this%Qtensor(2,2,vid) & 
            !                     & + this%PiBC(:,:,2,3)*this%Qtensor(2,3,vid) + this%PiBC(:,:,3,1)*this%Qtensor(3,1,vid) & 
            !                     & + this%PiBC(:,:,3,2)*this%Qtensor(3,2,vid) + this%PiBC(:,:,3,3)*this%Qtensor(3,3,vid))& 
            !                     & - (w(vid)*half*onebycsq)*(cx(vid)*this%Fx + cy(vid)*this%Fy + cz(vid)*this%Fz)
            
            this%fneqBC(:,:,vid) = (w(vid)*oneby2c4)*(this%PiBC(:,:,1,1)*this%Qtensor(1,1,vid) &
                                 & + two*this%PiBC(:,:,1,2)*this%Qtensor(1,2,vid) + two*this%PiBC(:,:,1,3)*this%Qtensor(1,3,vid) & 
                                 & + this%PiBC(:,:,2,2)*this%Qtensor(2,2,vid) + two*this%PiBC(:,:,2,3)*this%Qtensor(2,3,vid) & 
                                 & + this%PiBC(:,:,3,3)*this%Qtensor(3,3,vid))& 
                                 & - (w(vid)*half*onebycsq)*(cx(vid)*this%Fx + cy(vid)*this%Fy + cz(vid)*this%Fz)
        end do 


    end subroutine


    subroutine compute_f_BC(this, uBC, vBC, wBC, fBC)
        class(d3q19), intent(inout) :: this
        real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2)), intent(in) :: uBC, vBC, wBC
        real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),nvels), intent(in) :: fBC
        integer :: i, j, vid
        
        do vid = 1,nvels
            do j = 1,this%gp%zsz(2)
                !$omp simd
                do i = 1,this%gp%zsz(1)
                    call get_Feq_2ndOrder(uBC(i,j),vBC(i,j),wBC(i,j), this%rhoBC(i,j), &
                        & vid, this%Qtensor(:,:,vid),this%feqBC(i,j,vid))
                    
                    this%fneqBC(i,j,vid) = fBC(i,j,vid) - this%feqBC(i,j,vid)
                end do 
            end do
        end do 

    end subroutine 


    subroutine compute_rhoBC(this, fknown, wBC, fBC, sgForce)
        class(d3q19), intent(inout) :: this
        integer, dimension(:), intent(in) :: fknown
        real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2)), intent(in) :: wBC
        real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),nvels), intent(in) :: fBC
        integer, intent(in) :: sgForce
        integer :: vid
        
        this%rhoBC = zero
        do vid = 1,size(fparallel)
            this%rhoBC = this%rhoBC + fBC(:,:,fparallel(vid))
        end do 
        do vid = 1,size(fknown)
            this%rhoBC = this%rhoBC + two*fBC(:,:,fknown(vid))
        end do 
        this%rhoBC = (this%rhoBC + sgForce*half*this%Fz)/(one + wBC )

    end subroutine 
