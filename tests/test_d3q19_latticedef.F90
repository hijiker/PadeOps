#include "test_bolt_files/d3q19_hooks.F90"       

program test_d3q19_latticedef
    use d3q19mod, only: d3q19,test_lattice_definition
    use mpi
    use exits, only: message
    use kind_parameters, only: rkind

    implicit none

    type(d3q19) :: lattice 
    integer :: ierr
    real(rkind) :: sout

    call MPI_Init(ierr)               !<-- Begin MPI
    
    call test_lattice_definition(sout)

    if (abs(sout) < 1.d-14) then
        call message(0,"TEST PASSED.")
    else
        call message(0, "value of residue", sout)
        call message(0,"TEST FAILED.")
    end if 

    call MPI_Finalize(ierr)
end program 
