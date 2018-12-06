program test_MRT_matrixMult
    use kind_parameters, only: rkind
    use timer, only: tic, toc
    use MRT_D3Q19_coefficients, only: M, Minv
    use constants, only: zero 
    use random
    use mpi 
    implicit none

    integer, parameter :: nvels = 19
    integer, parameter :: nnodes = 1000000
    real(rkind), dimension(:,:), allocatable :: f, fnew
    real(rkind), dimension(nvels) :: ftmp, mvec
    integer :: i, ierr, iter
  
    call MPI_Init(ierr)

    allocate(f(nnodes,nvels))
    allocate(fnew(nnodes,nvels))

    call uniform_random(f,0.d0,1.d0,344)
    ! Do loops 
    call tic()
    do iter = 1,nnodes
        ftmp = f(iter,:)
        mvec = zero 
        do i = 1,nvels
            mvec = mvec+ ftmp(i)*M(:,i) 
        end do 

        ! Do collision

        ftmp = zero
        do i = 1,nvels
            ftmp = ftmp + mvec(i)*Minv(:,i) 
        end do
        fnew(iter,:) = ftmp
    end do 
    call toc()
    
    print*, "Error:", maxval(abs(fnew - f))

    call tic()
    do iter = 1,nnodes
        ftmp = f(iter,:)
        mvec = matmul(M,ftmp)

        ! Do collision
        ftmp = matmul(Minv,mvec)
        fnew(iter,:) = ftmp
    end do 
    call toc()
   
    print*, "Error:", maxval(abs(fnew - f))
    call MPI_Finalize(ierr)

end program 
