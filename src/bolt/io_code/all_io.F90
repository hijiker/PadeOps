subroutine read_inputfile(this, inputfile)
        class(d3q19), intent(inout) :: this
        character(len=*), intent(in) :: inputfile 

        integer :: nx, ny, nz, CollisionModel=0, restart_runID, restart_timeID, runID
        integer :: tid_vis, tid_restart, ierr, gradient_type = 1 
        character(len=clen) :: inputdir, outputdir
        logical ::  isZPeriodic=.false., useConstantBodyForce, useSGSmodel=.false.
        logical :: restartSimulation=.false., useSpaceTimeBodyForce = .false., restartWithTau = .false. 
        real(rkind) :: Re = 10.d0, delta_t = 0.1d0, delta_x = 0.2d0, Fx = zero, Fy = zero, Fz = zero, c_smag = 0.16d0

        namelist /INPUT/ nx, ny, nz, restartSimulation, restart_runID, restart_timeID, restartwithTau
        namelist /PHYSICS/ CollisionModel, useConstantBodyForce, useSGSmodel, isZperiodic, Re, delta_x, delta_t, Fx, Fy, Fz, useSpaceTimeBodyForce, gradient_type, c_smag 
        namelist /IO/ inputdir, outputdir, RunID, tid_vis, tid_restart 

        open(unit=123, file=trim(inputfile), form='FORMATTED', iostat=ierr)
        read(unit=123, NML=INPUT)
        read(unit=123, NML=PHYSICS)
        read(unit=123, NML=IO)
        close(123)

        this%CollisionModel = CollisionModel
        this%useSGSmodel = useSGSmodel
        this%isZPeriodic = isZPeriodic         
        this%inputdir = inputdir
        this%outputdir = outputdir
        this%useRestart = restartSimulation
        this%Re = Re
        this%useConstantBodyForce = useConstantBodyForce 
        this%restart_runID = restart_runID
        this%restart_timeID = restart_timeID
        this%useSpaceTimeBodyForce = useSpaceTimeBodyForce
        this%restartWithTau = restartwithTau
        this%gradient_type = gradient_type
        this%c_smag = c_smag

        this%nx = nx
        this%ny = ny
        this%nz = nz
     
        this%Fx = Fx
        this%Fy = Fy
        this%Fz = Fz

        this%delta_x = delta_x
        this%delta_t = delta_t

        this%RunID = RunID
        this%tid_vis = tid_vis
        this%tid_restart = tid_restart

end subroutine


subroutine dumpVisBuff(this,label)
    use decomp_2d_io
    class(d3q19), intent(in) :: this
    character(len=clen) :: tempname, fname
    character(len=4), intent(in) :: label

     write(tempname,"(A3,I2.2,A1,A4,A2,I6.6,A4)") "Run",this%runID, "_",label,"_t",this%step,".out"
     fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
     call decomp_2d_write_one(3,this%VisBuff,fname,this%gp)

end subroutine


subroutine dumpVisualizationFields(this)
    class(d3q19), intent(inout) :: this
    
    this%visBuff = this%rho   
    call this%dumpVisBuff("rhoF")
    
    this%visBuff = this%ux*this%delta_u   
    call this%dumpVisBuff("uVel")
    
    this%visBuff = this%uy*this%delta_u   
    call this%dumpVisBuff("vVel")
    
    this%visBuff = this%uz*this%delta_u   
    call this%dumpVisBuff("wVel")

    if (allocated(this%nusgs)) then
        this%visBuff = this%nusgs*this%delta_nu
        call this%dumpVisBuff("nSGS")
    end if 
end subroutine 

subroutine dumpRestart(this)
    use decomp_2d_io
    class(d3q19), intent(in) :: this
    character(len=clen) :: tempname, fname
    integer :: vid

    do vid = 1,nvels
        write(tempname,"(A7,A4,I2.2,A2,I3.3,A1,I6.6)") "RESTART", "_Run",this%runID, "_f",vid,".",this%step
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        call decomp_2d_write_one(3,this%f(:,:,:,vid),fname, this%gp)
    end do 

    write(tempname,"(A7,A4,I2.2,A4,A1,I6.6)") "RESTART", "_Run",this%runID, "_tau",".",this%step
    fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
    call decomp_2d_write_one(3,this%tau,fname, this%gp)
end subroutine


subroutine readRestart(this)
    use decomp_2d_io
    use mpi 
    use exits, only: gracefulExit
    class(d3q19), intent(inout) :: this
    character(len=clen) :: tempname, fname
    integer :: vid, ierr

    do vid = 1,nvels
        write(tempname,"(A7,A4,I2.2,A2,I3.3,A1,I6.6)") "RESTART", "_Run",this%restart_runID, "_f",vid,".",this%restart_timeID
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        open(777,file=trim(fname),status='old',iostat=ierr)
        if(ierr .ne. 0) then
            print*, "Rank:", nrank, ". File:", fname, " not found"
            call mpi_barrier(mpi_comm_world, ierr)
            call gracefulExit("File I/O issue.",44)
        end if 
        close(777)
        call decomp_2d_read_one(3,this%f(:,:,:,vid),fname, this%gp)
    end do 
    
    if (this%restartWithTau) then
        write(tempname,"(A7,A4,I2.2,A4,A1,I6.6)") "RESTART", "_Run",this%restart_runID, "_tau",".",this%restart_timeID
        fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        open(777,file=trim(fname),status='old',iostat=ierr)
        if(ierr .ne. 0) then
            print*, "Rank:", nrank, ". File:", fname, " not found"
            call mpi_barrier(mpi_comm_world, ierr)
            call gracefulExit("File I/O issue.",44)
        end if 
        call decomp_2d_read_one(3,this%tau,fname,this%gp)
    end if  

    this%step = this%restart_timeID

end subroutine 
