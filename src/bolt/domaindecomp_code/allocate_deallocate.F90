    subroutine allocate_lattice_memory(this)
        class(d3q19), intent(inout) :: this

        allocate(this%f  (this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3),nvels))
        allocate(this%ux (this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        allocate(this%uy (this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        allocate(this%uz (this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        allocate(this%rho(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        
        allocate(this%Fx (this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        allocate(this%Fy (this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        allocate(this%Fz (this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        
        allocate(this%tau(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        
        allocate(this%VisBuff(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))

        allocate(this%XYslice(this%gp%zsz(1),this%gp%zsz(2)))
        allocate(this%XZslice(this%gp%zsz(1),this%gp%zsz(3)))
        allocate(this%YZslice(this%gp%zsz(2),this%gp%zsz(3)))
        
        allocate(this%YZsendb(this%gp%zsz(2),this%gp%zsz(3)))
        allocate(this%XZsendb(this%gp%zsz(1),this%gp%zsz(3)))
        allocate(this%z1d_send(this%gp%zsz(3)))
        allocate(this%z1d_recv(this%gp%zsz(3)))

        if (.not. this%isZperiodic) then
            allocate(this%uTop(this%gp%zsz(1),this%gp%zsz(2)))
            allocate(this%vTop(this%gp%zsz(1),this%gp%zsz(2)))
            allocate(this%wTop(this%gp%zsz(1),this%gp%zsz(2)))
            
            allocate(this%uBot(this%gp%zsz(1),this%gp%zsz(2)))
            allocate(this%vBot(this%gp%zsz(1),this%gp%zsz(2)))
            allocate(this%wBot(this%gp%zsz(1),this%gp%zsz(2)))
            
            allocate(this%rhoBC(this%gp%zsz(1),this%gp%zsz(2)))
            allocate(this%PiBC(this%gp%zsz(1),this%gp%zsz(2),3,3))
        
            allocate(this%feqBC (this%gp%zsz(1),this%gp%zsz(2),nvels))
            allocate(this%fneqBC(this%gp%zsz(1),this%gp%zsz(2),nvels))
            
            allocate(this%tau_B(this%gp%zsz(1),this%gp%zsz(2)))
            allocate(this%tau_T(this%gp%zsz(1),this%gp%zsz(2)))
            
        end if 
   
        if (this%useSGSmodel) then
            allocate(this%nuSGS(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
            allocate(this%duidxj(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3),3,3))
            allocate(this%rbuffx1(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)))
            allocate(this%rbuffx2(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)))
            allocate(this%rbuffy1(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3)))
            allocate(this%rbuffy2(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3)))
        end if
        allocate(this%rbuffz1(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        allocate(this%rbuffz2(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))

        this%Fx = zero
        this%Fy = zero
        this%Fz = zero

    end subroutine 
    
    subroutine destroy(this)
        class(d3q19), intent(inout) :: this

        if (this%compute_stats) then
            call this%end_stats()
        end if 
        deallocate(this%f, this%ux, this%uy, this%uz, this%rho)
        deallocate(this%XYslice, this%XZslice, this%YZslice)
        deallocate(this%YZsendb, this%XZsendb)
        deallocate(this%z1d_send, this%z1d_recv)

    end subroutine 
