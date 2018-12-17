#include "test_bolt_files/d3q19_hooks.F90"       

program test_d3q19_uniform
    use d3q19mod, only: d3q19,test_lattice_definition
    use kind_parameters, only: rkind
    use mpi
    use reductions, only: p_maxval
    use exits, only: message
    use constants, only: one, ten
    use d3q19_testing, only: mvec_mrt_uniform, testing_type, utrue, vtrue, wtrue  
    implicit none

    type(d3q19) :: lattice 
    integer :: i, j, k, ierr, ii, jj, kk
    integer :: nx = 32, ny = 24, nz = 16

    real(rkind) :: tmp, accum

    call MPI_Init(ierr)               
    
    testing_type = 1

    lattice%nx = nx
    lattice%ny = ny
    lattice%nz = nz

    lattice%delta_x = one
    lattice%delta_t = one
    lattice%Re = ten 

    lattice%isZPeriodic = .true. 
    lattice%useConstantBodyForce = .false. 
    lattice%CollisionModel = 2

    call lattice%init(testing=.true.)

    call lattice%time_advance()

    accum = p_maxval(abs(lattice%ux*lattice%delta_u - utrue) &
          & + abs(lattice%uy*lattice%delta_u - vtrue) + abs(lattice%uz*lattice%delta_u - wtrue))

    call message(0,"ERROR ACCUMULATED:", accum)
    
    if (lattice%CollisionModel == 2) then
        accum = accum + p_maxval(maxval(abs(lattice%m_mrt - mvec_mrt_uniform)))
    end if 

    if (accum > 1.d-13) then
        call message(1,"TEST FAILED.")
    else
        call message(1,"TEST PASSED.")
    end if 

    call lattice%destroy()

    call MPI_Finalize(ierr)
end program 
