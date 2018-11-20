#include "channel3D_files/initialize.F90"       
#include "channel3D_files/temporalHook.F90"  

program channel3d
    use mpi
    use kind_parameters,  only: clen
    use temporalhook, only: doTemporalStuff
    use timer, only: tic, toc
    use d3q19mod, only: d3q19,test_lattice_definition

    type(d3q19) :: lattice 

    call test_lattice_definition()

end program 
