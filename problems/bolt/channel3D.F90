#include "channel3D_files/initialize.F90"       
#include "channel3D_files/temporalHook.F90"  

program channel3d
    use mpi
    use kind_parameters,  only: clen
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use d3q19mod, only: d3q19
    use bolt_hooks
    use d3q19_channel3D, only: BC_Type, Fbase 
    use get_initial_profiles_channel, only: utau, Retau 
    implicit none

    type(d3q19) :: lattice 
    character(len=clen) :: inputfile 
    integer :: ierr
    real(rkind) :: t_stop = 1.d-1 
    integer :: WallModelType = 0
    real(rkind) :: Retau_expect = 5.434960d02 
    real(rkind) :: utau_expect = 5.434960d-02 
    real(rkind) :: Forcex      = 2.953879d-03 
    namelist /CHANNEL_WM/ WallModelType, Retau_expect, utau_expect, Forcex

    call MPI_Init(ierr)               !<-- Begin MPI
    call GETARG(1,inputfile)          !<-- Get the location of the input file
    
    open(unit=123, file=trim(inputfile), form='FORMATTED', iostat=ierr)
    read(unit=123, NML=CHANNEL_WM)
    close(123)
    
    BC_Type = WallModelType
    Retau   = Retau_expect
    utau    = utau_expect
    Fbase   = Forcex

    call lattice%init(inputfile,.false.)

    call lattice%dumpVisualizationFields()

    call tic()
    do while (.not. lattice%EndSim)
        call lattice%time_advance()
        call doTemporalStuff(lattice)
    end do 
    call doTemporalStuff(lattice)
    call lattice%dumpVisualizationFields()

    call lattice%destroy()

    call MPI_Finalize(ierr)
end program 
