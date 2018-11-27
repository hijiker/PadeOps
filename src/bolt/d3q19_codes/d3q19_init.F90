    subroutine init(this,inputfile, testing)
        class(d3q19), intent(inout) :: this
        character(len=*), intent(in), optional :: inputfile 
        logical, optional, intent(in) :: testing 

        integer :: idx
        logical :: testing_mode

        if (present(testing)) then
            testing_mode = testing
        else
            testing_mode = .true. 
        end if 
        
        ! Populate the Qtensor
        do idx = 1,nvels
            this%Qtensor(1,1,idx) = cx(idx)*cx(idx) - csq
            this%Qtensor(2,2,idx) = cy(idx)*cy(idx) - csq
            this%Qtensor(3,3,idx) = cz(idx)*cz(idx) - csq

            this%Qtensor(1,2,idx) = cx(idx)*cy(idx) 
            this%Qtensor(2,1,idx) = cx(idx)*cy(idx) 

            this%Qtensor(1,3,idx) = cx(idx)*cz(idx) 
            this%Qtensor(3,1,idx) = cx(idx)*cz(idx) 

            this%Qtensor(2,3,idx) = cy(idx)*cz(idx) 
            this%Qtensor(3,2,idx) = cy(idx)*cz(idx) 
        end do 

        if (.not. testing_mode) then
            if (.not. present(inputfile)) then
                call GracefulExit("Input file was not passed in while initializing D3Q19",45)
            end if 

            call this%read_inputfile(inputfile)

        end if 
            
        call this%get_processor_topo()
        call this%allocate_lattice_memory()
        
        ! Compute lattice parameters 
        this%delta_u  = this%delta_x/this%delta_t
        this%delta_nu = (this%delta_x**2)/this%delta_t
        this%nu = (one/this%Re)/this%delta_nu
        this%tau = this%nu/csq + half 
        this%Fx = this%Fx*this%delta_t/this%delta_u
        this%Fy = this%Fy*this%delta_t/this%delta_u
        this%Fz = this%Fz*this%delta_t/this%delta_u


        if (this%useRestart) then
            call this%readRestart()   
            call this%compute_macroscopic()
        else
            call initfields_bolt(this%gp, inputfile, this%delta_x, this%rho, this%ux, this%uy, this%uz)
            this%ux = this%ux/this%delta_u 
            this%uy = this%uy/this%delta_u 
            this%uz = this%uz/this%delta_u 
            call this%initialize_f_to_feq()
            this%step = 0
        end if 

        call message(0, "D3Q19 lattice initialized.")
    end subroutine 


