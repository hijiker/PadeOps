subroutine read_inputfile(this, inputfile)
        class(d3q19), intent(inout) :: this
        character(len=*), intent(in) :: inputfile 

        integer :: nx, ny, nz, CollisionModel=0, restart_runID, restart_timeID, runID
        integer :: tid_vis, tid_restart, ierr, gradient_type = 1 
        character(len=clen) :: inputdir, outputdir, stats_dir
        logical ::  isZPeriodic=.false., useConstantBodyForce, useSGSmodel=.false., compute_stats = .false. 
        logical :: useBoundaryTau = .false., restartSimulation=.false., useSpaceTimeBodyForce = .false., restartWithTau = .false. 
        real(rkind) :: Re = 10.d0, delta_t = 0.1d0, delta_x = 0.2d0, Fx = zero, Fy = zero, Fz = zero, c_sgs = 0.16d0
        integer :: step_stop = 999999, tid_stats_dump, stats_freq, stats_start = 999999, sgs_model_type = 0

        namelist /INPUT/ nx, ny, nz, step_stop, restartSimulation, restart_runID, restart_timeID, restartwithTau
        namelist /PHYSICS/ CollisionModel, sgs_model_type,useConstantBodyForce, useSGSmodel, isZperiodic, Re, delta_x, delta_t, &
                    & Fx, Fy, Fz, useSpaceTimeBodyForce, gradient_type, c_sgs, useBoundaryTau  
        namelist /IO/ inputdir, outputdir, RunID, tid_vis, tid_restart 
        namelist /STATS/ compute_stats, tid_stats_dump, stats_freq, stats_dir,stats_start 


        open(unit=123, file=trim(inputfile), form='FORMATTED', iostat=ierr)
        read(unit=123, NML=INPUT)
        read(unit=123, NML=PHYSICS)
        read(unit=123, NML=IO)
        read(unit=123, NML=STATS)
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
        this%c_sgs = c_sgs
        this%sgs_model_type = sgs_model_type
        this%step_stop = step_stop
        this%useBoundaryTau = useBoundaryTau

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

        this%compute_stats = compute_stats
        this%tid_stats_dump = tid_stats_dump
        this%stats_freq = stats_freq
        this%stats_dir = stats_dir 
        this%stats_start = stats_start

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

subroutine dump_stats(this)
    use basic_io, only: write_2d_ascii
    use decomp_2d,only: nrank 
    class(d3q19), intent(inout) :: this
    character(len=clen) :: fname, tempname 

    write(tempname,"(A3,I2.2,A6,A2,I6.6,A2,I6.6,A4)") "Run",this%runid,"_stats","_t",this%step,"_n",this%dat_count,".stt"
    fname = this%stats_dir(:len_trim(this%stats_dir))//"/"//trim(tempname)
    
    this%dat_array = this%sum_array/real(this%dat_count,rkind)
    
    this%dat_array(:,4) = this%dat_array(:,4) - this%dat_array(:,1)*this%dat_array(:,1)
    this%dat_array(:,5) = this%dat_array(:,5) - this%dat_array(:,1)*this%dat_array(:,2)
    this%dat_array(:,6) = this%dat_array(:,6) - this%dat_array(:,1)*this%dat_array(:,3)
    this%dat_array(:,7) = this%dat_array(:,7) - this%dat_array(:,2)*this%dat_array(:,2)
    this%dat_array(:,8) = this%dat_array(:,8) - this%dat_array(:,2)*this%dat_array(:,3)
    this%dat_array(:,9) = this%dat_array(:,9) - this%dat_array(:,3)*this%dat_array(:,3)

    this%dat_array(:,1:3) = this%dat_array(:,1:3)*this%delta_u
    this%dat_array(:,4:9) = this%dat_array(:,4:9)*this%delta_u*this%delta_u
    
    ! Nothing to do to rho_mn (index: 10) and nusgs (index:11)
    
    if (nrank == 0) then
        call write_2d_ascii(this%dat_array, fname) 
    end if 
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

    if (this%useSGSmodel) then
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

    !write(tempname,"(A7,A4,I2.2,A4,A1,I6.6)") "RESTART", "_Run",this%runID, "_tau",".",this%step
    !fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
    !call decomp_2d_write_one(3,this%tau,fname, this%gp)
end subroutine


subroutine readRestart(this)
    use decomp_2d_io
    use mpi 
    use exits, only: message, gracefulExit
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
        !write(tempname,"(A7,A4,I2.2,A4,A1,I6.6)") "RESTART", "_Run",this%restart_runID, "_tau",".",this%restart_timeID
        !fname = this%OutputDir(:len_trim(this%OutputDir))//"/"//trim(tempname)
        !open(777,file=trim(fname),status='old',iostat=ierr)
        !if(ierr .ne. 0) then
        !    print*, "Rank:", nrank, ". File:", fname, " not found"
        !    call mpi_barrier(mpi_comm_world, ierr)
        !    call gracefulExit("File I/O issue.",44)
        !end if 
        !call decomp_2d_read_one(3,this%tau,fname,this%gp)
        call message(0,"WARNING: Restartwithtau redacted")
    end if  

    this%step = this%restart_timeID

end subroutine 
