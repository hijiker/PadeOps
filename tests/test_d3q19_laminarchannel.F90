#include "test_bolt_files/d3q19_hooks.F90"       

program test_d3q19_laminarchannel
    use d3q19mod, only: d3q19,test_lattice_definition
    use kind_parameters, only: rkind
    use mpi
    use reductions, only: p_maxval, p_minval
    use exits, only: message
    use constants, only: one, ten, two, three
    use d3q19_testing, only: testing_type, Reynolds_number, Force, utrue, vtrue, wtrue  
    use timer, only: tic, toc
    implicit none

    type(d3q19) :: lattice 
    integer :: ierr
    integer :: nx = 1, ny = 1, nz = 20
    real(rkind) :: Tstop = 150._rkind
    real(rkind) :: tmp, accum
    real(rkind) :: beta_t = 0.3_rkind

    real(rkind) :: umax_expect 

    call MPI_Init(ierr)               
    
    testing_type = 2

    lattice%nx = nx
    lattice%ny = ny
    lattice%nz = nz

    umax_expect = two*Force*Reynolds_number

    lattice%delta_x = two/real(nz-1,rkind)
    lattice%delta_t = beta_t*(lattice%delta_x/(sqrt(three)*umax_expect))
    lattice%Re = Reynolds_number
   

    lattice%isZPeriodic = .false. 
    lattice%useConstantBodyForce = .true.
    lattice%Fx = Force
    lattice%CollisionModel = 1

    call lattice%init(testing=.true.)

    call message(0, "Lattice details:")
    call message(1, "delta_x:",lattice%delta_x)
    call message(1, "delta_t:",lattice%delta_t)
    
    do while (lattice%GetPhysTime()<Tstop)
        call lattice%time_advance()
        if (mod(lattice%step,10000) == 0) then
            call message(0, "Time:", lattice%GetPhysTime())
            call message(1, "MaxU:", p_maxval(lattice%ux)*lattice%delta_u)
            call message(1, "MinU:", p_minval(lattice%ux)*lattice%delta_u)
            accum = p_maxval(abs(lattice%ux*lattice%delta_u - utrue) &
                    & + abs(lattice%uy*lattice%delta_u - vtrue) + abs(lattice%uz*lattice%delta_u - wtrue))
            call message(1,"Error:", accum)
            call message(1,"Step:", lattice%step)
        end if 
    end do 

    accum = p_maxval(abs(lattice%ux*lattice%delta_u - utrue) &
          & + abs(lattice%uy*lattice%delta_u - vtrue) + abs(lattice%uz*lattice%delta_u - wtrue))

    call lattice%compute_pi()

    call message("SUMMARY:")
    call message(1,"Time:",lattice%getPhysTime())
    call message(1,"Error:", accum)

    !print*, (one/lattice%delta_t)*lattice%Pitensor(3,1,1,:)/(lattice%tau(1,1,:)*(one/3.d0)*lattice%rho(1,1,:))

    if (accum > 1.7d-1) then
        call message(1,"TEST FAILED.")
    else
        call message(1,"TEST PASSED.")
    end if 

    call lattice%destroy()

    call MPI_Finalize(ierr)
end program 
