#include "test_bolt_files/d3q19_hooks.F90"       

program test_d3q19_laminarchannel
    use d3q19mod, only: d3q19,test_lattice_definition
    use kind_parameters, only: rkind
    use mpi
    use reductions, only: p_maxval
    use exits, only: message
    use constants, only: one, ten, two
    use d3q19_testing, only: testing_type, Reynolds_number, Force, utrue, vtrue, wtrue  
    use timer, only: tic, toc
    implicit none

    type(d3q19) :: lattice 
    integer :: ierr
    integer :: nx = 1, ny = 1, nz = 20
    real(rkind) :: Tstop = 50._rkind
    real(rkind) :: tmp, accum

    call MPI_Init(ierr)               
    
    testing_type = 2

    lattice%nx = nx
    lattice%ny = ny
    lattice%nz = nz

    lattice%delta_x = two/real(nz-1,rkind)
    lattice%delta_t = 0.25_rkind*lattice%delta_x
    lattice%Re = Reynolds_number
   
    lattice%isZPeriodic = .false. 
    lattice%useConstantBodyForce = .true.
    lattice%Fx = Force
    call lattice%init(testing=.true.)

    do while (lattice%GetPhysTime()<Tstop)
        call lattice%time_advance()
    end do 

    accum = p_maxval(abs(lattice%ux*lattice%delta_u - utrue) &
          & + abs(lattice%uy*lattice%delta_u - vtrue) + abs(lattice%uz*lattice%delta_u - wtrue))

    if (accum > 1.d-1) then
        call message(1,"TEST FAILED.")
    else
        call message(1,"TEST PASSED.")
    end if 

    call lattice%destroy()

    call MPI_Finalize(ierr)
end program 
