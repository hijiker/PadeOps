#include "test_bolt_files/d3q19_hooks.F90"       

program test_d3q19_streaming
    use d3q19mod, only: d3q19,velorder
    use kind_parameters, only: rkind
    use mpi
    use reductions, only: p_sum
    use exits, only: message
    use constants, only: one, ten

    implicit none

    type(d3q19) :: lattice 
    integer :: i, j, k, ierr, ii, jj, kk, vid
    integer :: nx = 32, ny = 24, nz = 16
    real(rkind), dimension(19) :: trueSum

    real(rkind) :: tmp, accum

    call MPI_Init(ierr)               !<-- Begin MPI
    
    lattice%nx = nx
    lattice%ny = ny
    lattice%nz = nz

    lattice%delta_x = one
    lattice%delta_t = one
    lattice%Re = ten 

    call lattice%init(testing=.true.)

    kk =  1
    do k = lattice%gp%zst(3),lattice%gp%zen(3)
        jj = 1
        do j = lattice%gp%zst(2),lattice%gp%zen(2)
            ii = 1
            do i = lattice%gp%zst(1),lattice%gp%zen(1)
                lattice%f(ii,jj,kk,:) = nx*ny*(k-1) + nx*(j-1) + i + 3.14d0 
                ii = ii  + 1
            end do 
            jj = jj + 1
        end do 
        kk = kk + 1
    end do 


    trueSum = [1.229393999589460d04, 1.228974294131363d04, 1.388285013669735d04, &
            &  1.243887559613250d04, 7.183356669572173d04, 1.773257145188483d04, &
            &  1.388879013259202d04, 1.244061853744609d04, 7.183950669161677d04, &
            &  1.773431439319851d04, 7.342841683241940d04, 1.788344704801725d04, &
            &  1.244481559202699d04, 1.388459307801110d04, 1.773851144777940d04, &
            &  7.183530963703620d04, 1.932742158858226d04, 7.198444229185481d04, &
            &  12288.d0]

    call lattice%stream()

    vid = 1
    accum = 0.d0 
    do while (vid < 20)
        tmp = p_sum(lattice%f(:,:,:,velorder(vid))/lattice%f(:,:,:,velorder(19)))
        call message(0, "Velocity index", vid)
        call message(1, "Error", abs(tmp - trueSum(vid))/tmp)
        call message("--------------------")
        accum = accum + abs(tmp - trueSum(vid))/tmp
        vid = vid + 1
    end do 

    if (accum > 1.d-13) then
        call message(0,"TEST FAILED.")
    else
        call message(0,"TEST PASSED.")
    end if 

    call lattice%destroy()

    call MPI_Finalize(ierr)
end program 
