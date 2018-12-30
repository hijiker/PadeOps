program Channel_Node2Cell
    use kind_parameters, only: rkind, clen
    use mpi 
    use decomp_2d
    use decomp_2d_io
    use exits, only: message, gracefulExit
    use PadeDerOps, only: Pade6stagg 
    implicit none
    
    character(len=clen) :: inputfile, outputdir, inputdir
    real(rkind), dimension(:,:,:), allocatable  :: fE, fC
    type(decomp_info) :: gpC, gpE
    type(Pade6stagg) :: der
    real(rkind) :: dz
    integer :: idx, nx, ny, nz, RID_input, RID_output, ierr, tid_input, tid_output, ioUnit
    character(len=clen) :: tempname, fname1, fname2
    integer, parameter :: nvels = 19
    
    namelist /INPUT/ nx, ny, nz, tid_input, tid_output, RID_input, RID_output, inputdir, outputdir
    
    
    call MPI_Init(ierr)               !<-- Begin MPI
    call GETARG(1,inputfile)          !<-- Get the location of the input file
    
    ioUnit = 11
    open(unit=ioUnit, file=trim(inputfile), form='FORMATTED')
    read(unit=ioUnit, NML=INPUT)
    close(ioUnit)
    
    
    call decomp_2d_init(nx, ny, nz, 0, 0)
    call get_decomp_info(gpC)
    call decomp_info_init(nx,ny,nz+1,gpE)
    dz = 2.d0/real(nz,rkind)
    
    allocate(fE(gpE%zsz(1), gpE%zsz(2), gpE%zsz(3)))
    allocate(fC(gpC%zsz(1), gpC%zsz(2), gpC%zsz(3)))
    
    do idx = 1,nvels
         write(tempname,"(A7,A4,I2.2,A2,I3.3,A1,I6.6)") "RESTART", "_Run",RID_input, "_f",idx,".",tid_input
         fname1 = InputDir(:len_trim(InputDir))//"/"//trim(tempname)
         open(777,file=trim(fname1),status='old',iostat=ierr)
         if(ierr .ne. 0) then
             print*, "Rank:", nrank, ". File:", fname1, " not found"
             call mpi_barrier(mpi_comm_world, ierr)
             call gracefulExit("File I/O issue.",44)
         end if 
         close(777)
         call decomp_2d_read_one(3,fE,fname1, gpE)
    
         call der%init(gpC, gpC, gpE, gpE, dz, 1, .false.)
         call der%interpz_E2C(fE, fC, 0, 0)
          
         write(tempname,"(A7,A4,I2.2,A2,I3.3,A1,I6.6)") "RESTART", "_Run",RID_output, "_f",idx,".",tid_output
         fname1 = OutputDir(:len_trim(OutputDir))//"/"//trim(tempname)
         call decomp_2d_write_one(3,fC,fname2, gpC)

         call message(0, "Restart file generated for velocity index", idx)
    end do 
    
   call message(0, "Conversion complete.")
   call MPI_Finalize(ierr)   

end program 
