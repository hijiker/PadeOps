#include "channelDNS_files/initialize.F90"       
#include "channelDNS_files/temporalHook.F90"  

program channelDNS
    use mpi
    use kind_parameters,  only: clen
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use d3q19mod, only: d3q19
    use bolt_hooks

    implicit none

    type(d3q19) :: lattice 
    character(len=clen) :: inputfile 
    integer :: ierr

    call MPI_Init(ierr)               !<-- Begin MPI
    call GETARG(1,inputfile)          !<-- Get the location of the input file

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
