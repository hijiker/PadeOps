    subroutine stream(this)
        use mpi 
        use kind_parameters, only: mpirkind
        class(d3q19), intent(inout) :: this
        integer :: vid, nx, ny, nz, ierr
        integer :: send_req1, recv_req1, tag=123, status(MPI_STATUS_SIZE)
        integer :: send_req2, recv_req2

        nx = this%gp%zsz(1)
        ny = this%gp%zsz(2)
        nz = this%gp%zsz(3)

        ! WARNING: This subroutine has not been optimized for cache usage.
        
        ! Population 1:
        vid = velOrder(1)
        call MPI_IRECV(this%YZslice,ny*nz,mpirkind,this%XneighLeft,tag,MPI_COMM_WORLD,recv_req1,ierr)
        this%YZsendb = this%f(nx,:,:,vid)
        call MPI_ISEND(this%YZsendb,ny*nz,mpirkind,this%XneighRight,tag,MPI_COMM_WORLD,send_req1,ierr)
        this%f(2:nx,:,:,vid) = this%f(1:nx-1,:,:,vid)
        call MPI_WAIT(recv_req1,status,ierr)
        this%f(1,:,:,vid) = this%YZslice
        call MPI_WAIT(send_req1,status,ierr)
        call mpi_barrier(mpi_comm_world, ierr)

        ! Population 2: 
        vid = velOrder(2)
        call MPI_IRECV(this%YZslice,ny*nz,mpirkind,this%XneighRight,tag,MPI_COMM_WORLD,recv_req1,ierr)
        this%YZsendb = this%f(1,:,:,vid)
        call MPI_ISEND(this%YZsendb,ny*nz,mpirkind,this%XneighLeft,tag,MPI_COMM_WORLD,send_req1,ierr)
        this%f(1:nx-1,:,:,vid) = this%f(2:nx,:,:,vid)
        call MPI_WAIT(recv_req1,status,ierr)
        this%f(nx,:,:,vid) = this%YZslice
        call MPI_WAIT(send_req1,status,ierr)
        call mpi_barrier(mpi_comm_world, ierr)

        ! Population 3: 
        vid = velOrder(3)
        call MPI_IRECV(this%XZslice,nx*nz,mpirkind,this%YneighDown,tag,MPI_COMM_WORLD,recv_req1,ierr)
        this%XZsendb = this%f(:,ny,:,vid)
        call MPI_ISEND(this%XZsendb,nx*nz,mpirkind,this%YneighUp,tag,MPI_COMM_WORLD,send_req1,ierr)
        this%f(:,2:ny,:,vid) = this%f(:,1:ny-1,:,vid)
        call MPI_WAIT(recv_req1,status,ierr)
        this%f(:,1,:,vid) = this%XZslice
        call MPI_WAIT(send_req1,status,ierr)
        call mpi_barrier(mpi_comm_world, ierr)

        ! Population 4: 
        vid = velOrder(4)
        call MPI_IRECV(this%XZslice,nx*nz,mpirkind,this%YneighUp,tag,MPI_COMM_WORLD,recv_req1,ierr)
        this%XZsendb = this%f(:,1,:,vid)
        call MPI_ISEND(this%XZsendb,nx*nz,mpirkind,this%YneighDown,tag,MPI_COMM_WORLD,send_req1,ierr)
        this%f(:,1:ny-1,:,vid) = this%f(:,2:ny,:,vid)
        call MPI_WAIT(recv_req1,status,ierr)
        this%f(:,ny,:,vid) = this%XZslice
        call MPI_WAIT(send_req1,status,ierr)
        call mpi_barrier(mpi_comm_world, ierr)

        ! Population 5: 
        vid = velOrder(5)
        this%XYslice = this%f(:,:,nz,vid)
        this%f(:,:,2:nz,vid) = this%f(:,:,1:nz-1,vid)
        this%f(:,:,1,vid) = this%XYslice

        ! Population 6: 
        vid = velOrder(6)
        this%XYslice = this%f(:,:,1,vid)
        this%f(:,:,1:nz-1,vid) = this%f(:,:,2:nz,vid)
        this%f(:,:,nz,vid) = this%XYslice

        ! Population 7: 
        vid = velOrder(7)
        call MPI_IRECV(this%YZslice,ny*nz,mpirkind,this%XneighLeft,tag,MPI_COMM_WORLD,recv_req1,ierr)
        this%YZsendb = this%f(nx,:,:,vid)
        call MPI_ISEND(this%YZsendb,ny*nz,mpirkind,this%XneighRight,tag,MPI_COMM_WORLD,send_req1,ierr)
        call MPI_IRECV(this%XZslice,nx*nz,mpirkind,this%YneighDown,tag,MPI_COMM_WORLD,recv_req2,ierr)
        this%XZsendb = this%f(:,ny,:,vid)
        call MPI_ISEND(this%XZsendb,nx*nz,mpirkind,this%YneighUp,tag,MPI_COMM_WORLD,send_req2,ierr)
        this%f(2:nx,2:ny,:,vid) = this%f(1:nx-1,1:ny-1,:,vid)
        call MPI_WAIT(recv_req1,status,ierr)
        this%f(1,2:ny,:,vid) = this%YZslice(1:ny-1,:)
        call MPI_WAIT(recv_req2,status,ierr)
        this%f(2:nx,1,:,vid) = this%XZslice(1:nx-1,:)
        call MPI_WAIT(send_req1,status,ierr)
        call MPI_WAIT(send_req2,status,ierr)
        call MPI_IRECV(this%z1d_recv,nz,mpirkind,this%XneighLeft,tag,MPI_COMM_WORLD,recv_req1,ierr)
        this%z1d_send = this%XZslice(nx,:)
        call MPI_ISEND(this%z1d_send,nz,mpirkind,this%XneighRight,tag,MPI_COMM_WORLD,send_req1,ierr)
        call MPI_WAIT(recv_req1,status,ierr)
        this%f(1,1,:,vid) = this%z1d_recv 
        call MPI_WAIT(send_req1,status,ierr)
        call mpi_barrier(mpi_comm_world, ierr)

        ! Population 8: 
        vid = velOrder(8)
        call MPI_IRECV(this%YZslice,ny*nz,mpirkind,this%XneighRight,tag,MPI_COMM_WORLD,recv_req1,ierr)
        this%YZsendb = this%f(1,:,:,vid)
        call MPI_ISEND(this%YZsendb,ny*nz,mpirkind,this%XneighLeft,tag,MPI_COMM_WORLD,send_req1,ierr)
        call MPI_IRECV(this%XZslice,nx*nz,mpirkind,this%YneighUp,tag,MPI_COMM_WORLD,recv_req2,ierr)
        this%XZsendb = this%f(:,1,:,vid)
        call MPI_ISEND(this%XZsendb,nx*nz,mpirkind,this%YneighDown,tag,MPI_COMM_WORLD,send_req2,ierr)
        this%f(1:nx-1,1:ny-1,:,vid) = this%f(2:nx,2:ny,:,vid)
        call MPI_WAIT(recv_req1,status,ierr)
        this%f(nx,1:ny-1,:,vid) = this%YZslice(2:ny,:)
        call MPI_WAIT(recv_req2,status,ierr)
        this%f(1:nx-1,ny,:,vid) = this%XZslice(2:nx,:)
        call MPI_WAIT(send_req1,status,ierr)
        call MPI_WAIT(send_req2,status,ierr)
        call MPI_IRECV(this%z1d_recv,nz,mpirkind,this%XneighRight,tag,MPI_COMM_WORLD,recv_req1,ierr)
        this%z1d_send = this%XZslice(1,:)
        call MPI_ISEND(this%z1d_send,nz,mpirkind,this%XneighLeft,tag,MPI_COMM_WORLD,send_req1,ierr)
        call MPI_WAIT(recv_req1,status,ierr)
        this%f(nx,ny,:,vid) = this%z1d_recv 
        call MPI_WAIT(send_req1,status,ierr)
        call mpi_barrier(mpi_comm_world, ierr)

        ! Population 9: 
        vid = velOrder(9)
        call MPI_IRECV(this%YZslice,ny*nz,mpirkind,this%XneighLeft,tag,MPI_COMM_WORLD,recv_req1,ierr)
        this%YZsendb = this%f(nx,:,:,vid)
        call MPI_ISEND(this%YZsendb,ny*nz,mpirkind,this%XneighRight,tag,MPI_COMM_WORLD,send_req1,ierr)
        this%XYslice = this%f(:,:,nz,vid)
        this%f(2:nx,:,2:nz,vid) = this%f(1:nx-1,:,1:nz-1,vid)
        this%f(2:nx,:,1,vid) = this%XYslice(1:nx-1,:)
        call MPI_WAIT(recv_req1,status,ierr)
        this%f(1,:,2:nz,vid) = this%YZslice(:,1:nz-1)
        this%f(1,:,1,vid) = this%YZslice(:,nz)
        call MPI_WAIT(send_req1,status,ierr)
        call mpi_barrier(mpi_comm_world, ierr)

        ! Population 10: 
        vid = velOrder(10)
        call MPI_IRECV(this%YZslice,ny*nz,mpirkind,this%XneighRight,tag,MPI_COMM_WORLD,recv_req1,ierr)
        this%YZsendb = this%f(1,:,:,vid)
        call MPI_ISEND(this%YZsendb,ny*nz,mpirkind,this%XneighLeft,tag,MPI_COMM_WORLD,send_req1,ierr)
        this%XYslice = this%f(:,:,1,vid)
        this%f(1:nx-1,:,1:nz-1,vid) = this%f(2:nx,:,2:nz,vid)
        this%f(1:nx-1,:,nz,vid) = this%XYslice(2:nx,:)
        call MPI_WAIT(recv_req1,status,ierr)
        this%f(nx,:,1:nz-1,vid) = this%YZslice(:,2:nz)
        this%f(nx,:,nz,vid) = this%YZslice(:,1)
        call MPI_WAIT(send_req1,status,ierr)
        call mpi_barrier(mpi_comm_world, ierr)

        ! Population 11: 
        vid = velOrder(11)
        call MPI_IRECV(this%XZslice,nx*nz,mpirkind,this%YneighDown,tag,MPI_COMM_WORLD,recv_req1,ierr)
        this%XZsendb = this%f(:,ny,:,vid)
        call MPI_ISEND(this%XZsendb,nx*nz,mpirkind,this%YneighUp,tag,MPI_COMM_WORLD,send_req1,ierr)
        this%XYslice = this%f(:,:,nz,vid)
        this%f(:,2:ny,2:nz,vid) = this%f(:,1:ny-1,1:nz-1,vid)
        this%f(:,2:ny,1,vid) = this%XYslice(:,1:ny-1)
        call MPI_WAIT(recv_req1,status,ierr)
        this%f(:,1,2:nz,vid) = this%XZslice(:,1:nz-1)
        this%f(:,1,1,vid) = this%XZslice(:,nz)
        call MPI_WAIT(send_req1,status,ierr)
        call mpi_barrier(mpi_comm_world, ierr)

        ! Population 12: 
        vid = velOrder(12)
        call MPI_IRECV(this%XZslice,nx*nz,mpirkind,this%YneighUp,tag,MPI_COMM_WORLD,recv_req1,ierr)
        this%XZsendb = this%f(:,1,:,vid)
        call MPI_ISEND(this%XZsendb,nx*nz,mpirkind,this%YneighDown,tag,MPI_COMM_WORLD,send_req1,ierr)
        this%XYslice = this%f(:,:,1,vid)
        this%f(:,1:ny-1,1:nz-1,vid) = this%f(:,2:ny,2:nz,vid)
        this%f(:,1:ny-1,nz,vid) = this%XYslice(:,2:ny)
        call MPI_WAIT(recv_req1,status,ierr)
        this%f(:,ny,1:nz-1,vid) = this%XZslice(:,2:nz)
        this%f(:,ny,nz,vid) = this%XZslice(:,1)
        call MPI_WAIT(send_req1,status,ierr)
        call mpi_barrier(mpi_comm_world, ierr)

        ! Population 13: 
        vid = velOrder(13)
        call MPI_IRECV(this%YZslice,ny*nz,mpirkind,this%XneighLeft,tag,MPI_COMM_WORLD,recv_req1,ierr)
        this%YZsendb = this%f(nx,:,:,vid)
        call MPI_ISEND(this%YZsendb,ny*nz,mpirkind,this%XneighRight,tag,MPI_COMM_WORLD,send_req1,ierr)
        call MPI_IRECV(this%XZslice,nx*nz,mpirkind,this%YneighUp,tag,MPI_COMM_WORLD,recv_req2,ierr)
        this%XZsendb = this%f(:,1,:,vid)
        call MPI_ISEND(this%XZsendb,nx*nz,mpirkind,this%YneighDown,tag,MPI_COMM_WORLD,send_req2,ierr)
        this%f(2:nx,1:ny-1,:,vid) = this%f(1:nx-1,2:ny,:,vid)
        call MPI_WAIT(recv_req1,status,ierr)
        this%f(1,1:ny-1,:,vid) = this%YZslice(2:ny,:)
        call MPI_WAIT(recv_req2,status,ierr)
        this%f(2:nx,ny,:,vid) = this%XZslice(1:nx-1,:)
        call MPI_WAIT(send_req1,status,ierr)
        call MPI_WAIT(send_req2,status,ierr)
        call MPI_IRECV(this%z1d_recv,nz,mpirkind,this%XneighLeft,tag,MPI_COMM_WORLD,recv_req1,ierr)
        this%z1d_send = this%XZslice(nx,:)
        call MPI_ISEND(this%z1d_send,nz,mpirkind,this%XneighRight,tag,MPI_COMM_WORLD,send_req1,ierr)
        call MPI_WAIT(recv_req1,status,ierr)
        this%f(1,ny,:,vid) = this%z1d_recv 
        call MPI_WAIT(send_req1,status,ierr)
        call mpi_barrier(mpi_comm_world, ierr)
        

        ! Population 14: 
        vid = velOrder(14)
        call MPI_IRECV(this%YZslice,ny*nz,mpirkind,this%XneighRight,tag,MPI_COMM_WORLD,recv_req1,ierr)
        this%YZsendb = this%f(1,:,:,vid)
        call MPI_ISEND(this%YZsendb,ny*nz,mpirkind,this%XneighLeft,tag,MPI_COMM_WORLD,send_req1,ierr)
        call MPI_IRECV(this%XZslice,nx*nz,mpirkind,this%YneighDown,tag,MPI_COMM_WORLD,recv_req2,ierr)
        this%XZsendb = this%f(:,ny,:,vid)
        call MPI_ISEND(this%XZsendb,nx*nz,mpirkind,this%YneighUp,tag,MPI_COMM_WORLD,send_req2,ierr)
        this%f(1:nx-1,2:ny,:,vid) = this%f(2:nx,1:ny-1,:,vid)
        call MPI_WAIT(recv_req1,status,ierr)
        this%f(nx,2:ny,:,vid) = this%YZslice(1:ny-1,:)
        call MPI_WAIT(recv_req2,status,ierr)
        this%f(1:nx-1,1,:,vid) = this%XZslice(2:nx,:)
        call MPI_IRECV(this%z1d_recv,nz,mpirkind,this%XneighRight,tag,MPI_COMM_WORLD,recv_req1,ierr)
        this%z1d_send = this%XZslice(1,:)
        call MPI_WAIT(send_req1,status,ierr)
        call MPI_ISEND(this%z1d_send,nz,mpirkind,this%XneighLeft,tag,MPI_COMM_WORLD,send_req1,ierr)
        call MPI_WAIT(recv_req1,status,ierr)
        this%f(nx,1,:,vid) = this%z1d_recv 
        call MPI_WAIT(send_req1,status,ierr)
        call MPI_WAIT(send_req2,status,ierr)
        call mpi_barrier(mpi_comm_world, ierr)
        

        ! Population 15:
        vid = velOrder(15)
        call MPI_IRECV(this%YZslice,ny*nz,mpirkind,this%XneighLeft,tag,MPI_COMM_WORLD,recv_req1,ierr)
        this%YZsendb = this%f(nx,:,:,vid)
        call MPI_ISEND(this%YZsendb,ny*nz,mpirkind,this%XneighRight,tag,MPI_COMM_WORLD,send_req1,ierr)
        this%XYslice = this%f(:,:,1,vid)
        this%f(2:nx,:,1:nz-1,vid) = this%f(1:nx-1,:,2:nz,vid)
        this%f(2:nx,:,nz,vid) = this%XYslice(1:nx-1,:)
        call MPI_WAIT(recv_req1,status,ierr)
        this%f(1,:,1:nz-1,vid) = this%YZslice(:,2:nz) 
        this%f(1,:,nz,vid) = this%YZslice(:,1)
        call MPI_WAIT(send_req1,status,ierr)
        call mpi_barrier(mpi_comm_world, ierr)

        ! Population 16:
        vid = velOrder(16)
        call MPI_IRECV(this%YZslice,ny*nz,mpirkind,this%XneighRight,tag,MPI_COMM_WORLD,recv_req1,ierr)
        this%YZsendb = this%f(1,:,:,vid)
        call MPI_ISEND(this%YZsendb,ny*nz,mpirkind,this%XneighLeft,tag,MPI_COMM_WORLD,send_req1,ierr)
        this%XYslice = this%f(:,:,nz,vid)
        this%f(1:nx-1,:,2:nz,vid) = this%f(2:nx,:,1:nz-1,vid)
        this%f(1:nx-1,:,1,vid) = this%XYslice(2:nx,:)
        call MPI_WAIT(recv_req1,status,ierr)
        this%f(nx,:,2:nz,vid) = this%YZslice(:,1:nz-1) 
        this%f(nx,:,1,vid) = this%YZslice(:,nz)
        call MPI_WAIT(send_req1,status,ierr)
        call mpi_barrier(mpi_comm_world, ierr)

        ! Population 17:
        vid = velOrder(17)
        call MPI_IRECV(this%XZslice,nx*nz,mpirkind,this%YneighDown,tag,MPI_COMM_WORLD,recv_req1,ierr)
        this%XZsendb = this%f(:,ny,:,vid)
        call MPI_ISEND(this%XZsendb,nx*nz,mpirkind,this%YneighUp,tag,MPI_COMM_WORLD,send_req1,ierr)
        this%XYslice = this%f(:,:,1,vid)
        this%f(:,2:ny,1:nz-1,vid) = this%f(:,1:ny-1,2:nz,vid)
        this%f(:,2:ny,nz,vid) = this%XYslice(:,1:ny-1)
        call MPI_WAIT(recv_req1,status,ierr)
        this%f(:,1,1:nz-1,vid) = this%XZslice(:,2:nz)
        this%f(:,1,nz,vid) = this%XZslice(:,1)
        call MPI_WAIT(send_req1,status,ierr)
        call mpi_barrier(mpi_comm_world, ierr)

        ! Population 18:
        vid = velOrder(18)
        call MPI_IRECV(this%XZslice,nx*nz,mpirkind,this%YneighUp,tag,MPI_COMM_WORLD,recv_req1,ierr)
        this%XZsendb = this%f(:,1,:,vid)
        call MPI_ISEND(this%XZsendb,nx*nz,mpirkind,this%YneighDown,tag,MPI_COMM_WORLD,send_req1,ierr)
        this%XYslice = this%f(:,:,nz,vid)
        this%f(:,1:ny-1,2:nz,vid) = this%f(:,2:ny,1:nz-1,vid)
        this%f(:,1:ny-1,1,vid) = this%XYslice(:,2:ny)
        call MPI_WAIT(recv_req1,status,ierr)
        this%f(:,ny,2:nz,vid) = this%XZslice(:,1:nz-1)
        this%f(:,ny,1,vid) = this%XZslice(:,nz)
        call MPI_WAIT(send_req1,status,ierr)
        call mpi_barrier(mpi_comm_world, ierr)

        ! Population 19:
        ! Nothing to do here. 

    end subroutine 
