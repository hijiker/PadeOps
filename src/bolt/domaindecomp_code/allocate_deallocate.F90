    subroutine allocate_lattice_memory(this)
        class(d3q19), intent(inout) :: this

        allocate(this%f  (this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3),nvels))
        allocate(this%ux (this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        allocate(this%uy (this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        allocate(this%uz (this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        allocate(this%rho(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        
        allocate(this%tau(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        allocate(this%nu(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        
        allocate(this%VisBuff(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))

        allocate(this%XYslice(this%gp%zsz(1),this%gp%zsz(2)))
        allocate(this%XZslice(this%gp%zsz(1),this%gp%zsz(3)))
        allocate(this%YZslice(this%gp%zsz(2),this%gp%zsz(3)))
        
        allocate(this%YZsendb(this%gp%zsz(2),this%gp%zsz(3)))
        allocate(this%XZsendb(this%gp%zsz(1),this%gp%zsz(3)))
        allocate(this%z1d_send(this%gp%zsz(3)))
        allocate(this%z1d_recv(this%gp%zsz(3)))

        this%Fx = zero
        this%Fy = zero
        this%Fz = zero

        if (.not. this%isZperiodic) then
            allocate(this%uTop(this%gp%zsz(1),this%gp%zsz(2)))
            allocate(this%vTop(this%gp%zsz(1),this%gp%zsz(2)))
            allocate(this%wTop(this%gp%zsz(1),this%gp%zsz(2)))
            
            allocate(this%uBot(this%gp%zsz(1),this%gp%zsz(2)))
            allocate(this%vBot(this%gp%zsz(1),this%gp%zsz(2)))
            allocate(this%wBot(this%gp%zsz(1),this%gp%zsz(2)))
            
            allocate(this%rhoBC(this%gp%zsz(1),this%gp%zsz(2)))
            allocate(this%PiBC(this%gp%zsz(1),this%gp%zsz(2),6))
        
            allocate(this%feqBC (this%gp%zsz(1),this%gp%zsz(2),9))
            allocate(this%fneqBC(this%gp%zsz(1),this%gp%zsz(2),9))
            
        end if 
    end subroutine 
    
    subroutine destroy(this)
        class(d3q19), intent(inout) :: this

        deallocate(this%f, this%ux, this%uy, this%uz, this%rho)
        deallocate(this%XYslice, this%XZslice, this%YZslice)
        deallocate(this%YZsendb, this%XZsendb)
        deallocate(this%z1d_send, this%z1d_recv)

    end subroutine 
