#include "channel3D_files/initialize.F90"       
#include "channel3D_files/temporalHook.F90"  

program channel3d
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
    real(rkind) :: t_stop = 1.d-1 

    call MPI_Init(ierr)               !<-- Begin MPI
    call GETARG(1,inputfile)          !<-- Get the location of the input file

    call lattice%init(inputfile,.false.)

    call lattice%dumpVisualizationFields()

    call tic()
    do while (lattice%step<100)
        call lattice%time_advance()
        call doTemporalStuff(lattice)
    end do 
    call doTemporalStuff(lattice)
    call lattice%dumpVisualizationFields()

    call lattice%destroy()

    call MPI_Finalize(ierr)
end program 
