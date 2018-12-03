#include "test_bolt_files/d3q19_hooks.F90"       

program test_d3q19_womersley
    use d3q19mod, only: d3q19,test_lattice_definition
    use kind_parameters, only: rkind
    use mpi
    use reductions, only: p_maxval, p_minval
    use exits, only: message
    use constants, only: one, ten, two, three
    use d3q19_testing, only: testing_type, Reynolds_number, p0, omega, utrue, vtrue, wtrue, update_utrue_wormersley 
    use timer, only: tic, toc
    implicit none

    type(d3q19) :: lattice 
    integer :: ierr
    integer :: nx = 1, ny = 1, nz = 80
    real(rkind) :: Tstop = 10._rkind
    real(rkind) :: tmp, accum
    real(rkind) :: beta_t = 0.2_rkind
    real(rkind) :: ForceX
    real(rkind) :: umax_expect 
    integer :: stop_step = 1250

    call MPI_Init(ierr)               
    
    testing_type = 3

    lattice%nx = nx
    lattice%ny = ny
    lattice%nz = nz

    Reynolds_number = one/(0.3472237895076d0)
    umax_expect = p0*Reynolds_number/(4._rkind) 
    

    lattice%delta_x = two/real(nz-1,rkind)
    lattice%delta_t = beta_t*(lattice%delta_x/(sqrt(three)*umax_expect))
    lattice%Re = Reynolds_number
   

    lattice%CollisionModel = 1
    lattice%isZPeriodic = .false. 
    lattice%useConstantBodyForce = .true.
    lattice%Fx = p0
    call lattice%init(testing=.true.)


    call message(0, "Lattice details:")
    call message(1, "delta_x:",lattice%delta_x)
    call message(1, "delta_t:",lattice%delta_t)
   
    !do while (lattice%GetPhysTime()<Tstop)
    do while (lattice%step<stop_step)
        ForceX = p0*cos(omega*lattice%GetPhysTime())
        call lattice%update_bodyForce(ForceX)
        call lattice%time_advance()
        if (mod(lattice%step,10000) == 0) then
            call update_utrue_wormersley(lattice%GetPhysTime())
            call message(0, "Time:", lattice%GetPhysTime())
            call message(1, "MaxU:", p_maxval(lattice%ux)*lattice%delta_u)
            call message(1, "MinU:", p_minval(lattice%ux)*lattice%delta_u)
            accum = p_maxval(abs(lattice%ux*lattice%delta_u - utrue) &
                    & + abs(lattice%uy*lattice%delta_u - vtrue) + abs(lattice%uz*lattice%delta_u - wtrue))
            call message(1,"Error:", accum)
            call message(1,"Step:", lattice%step)
        end if 
    end do 

    call update_utrue_wormersley(lattice%GetPhysTime())
    accum = p_maxval(abs(lattice%ux*lattice%delta_u - utrue) &
          & + abs(lattice%uy*lattice%delta_u - vtrue) + abs(lattice%uz*lattice%delta_u - wtrue))

   
    call message("SUMMARY:")
    call message(1,"Time:",lattice%GetPhysTime()) 
    call message(1,"Steps:",lattice%step) 
    call message(1,"Error:",accum) 
    if (accum > 1.d-1) then
        call message(1,"TEST FAILED.")
    else
        call message(1,"TEST PASSED.")
    end if 

    call lattice%destroy()

    call MPI_Finalize(ierr)
end program 
