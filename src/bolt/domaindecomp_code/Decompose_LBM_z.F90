    subroutine get_processor_topo(this)
        use mpi
        class(d3q19), intent(inout) :: this
        logical, dimension(3) :: periodicbcs
        integer :: nrow, ncol, ntasks, ierr 
        integer :: DECOMP_CART_Z
    
        ! Get processor decomposition
        call MPI_COMM_SIZE(MPI_COMM_WORLD,ntasks,ierr) 
        call get_proc_factors(ntasks,nrow,ncol)
        if (nrow*ncol .ne. ntasks) then
            call GracefulExit("Processor decomposition failed.", 123)
        end if 

        ! Set all BCs to true (will handle non-periodicity separately)
        periodicbcs(1) = .true.; periodicbcs(2) = .true.; periodicbcs(3) = .true.
        call decomp_2d_init(this%nx, this%ny, this%nz, nrow, ncol, periodicbcs)
        call get_decomp_info(this%gp)

        ! Get neighbors
        call MPI_CART_CREATE(MPI_COMM_WORLD,2,[nrow,ncol],[.true.,.true.],.false.,DECOMP_CART_Z,ierr)
        call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Z, 0, 1, this%XneighLeft, this%XneighRight, ierr) ! east & west
        call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Z, 1, 1, this%YneighDown, this%YneighUp, ierr) ! east & west

        call message(0,"Processor decomposition information:")
        call message(1,"Nrows:", nrow)
        call message(1,"Ncols:", ncol)

    end subroutine 


    pure subroutine get_proc_factors(nproc, nrow, ncol)
        integer, intent(in) :: nproc
        integer, intent(out) :: ncol, nrow

        integer :: res

        nrow = floor(sqrt(real(nproc)))
        res = mod(nproc,nrow)

        do while ((res .ne. 0) .and. (nrow>1)) 
            nrow = nrow - 1
            res = mod(nproc,nrow)
        end do

        ncol = nproc/nrow

    end subroutine 
