subroutine compute_duidxj(this, dimensionalize)
    class(d3q19), intent(inout) :: this
    logical, intent(in) :: dimensionalize
    real(rkind) :: oneBydeltaT


    call this%get_gradient(this%ux,this%duidxj(:,:,:,1,1),this%duidxj(:,:,:,1,2),this%duidxj(:,:,:,1,3))
    call this%get_gradient(this%uy,this%duidxj(:,:,:,2,1),this%duidxj(:,:,:,2,2),this%duidxj(:,:,:,2,3))
    call this%get_gradient(this%uz,this%duidxj(:,:,:,3,1),this%duidxj(:,:,:,3,2),this%duidxj(:,:,:,3,3))
   
    if (dimensionalize) then
        oneByDeltaT = 1.d0/this%delta_t

        this%duidxj = this%duidxj*oneByDeltaT
    end if 
end subroutine

subroutine get_gradient(this, u, dudx, dudy, dudz)
    class(d3q19), intent(inout) :: this
    real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)), intent(in ) :: u
    real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)), intent(out) :: dudx, dudy, dudz

    call this%get_ddz(u,dudz)

    call transpose_z_to_y(u,this%rbuffy1,this%gp)
    call this%get_ddy(this%rbuffy1,this%rbuffy2)
    call transpose_y_to_z(this%rbuffy2,dudy)

    call transpose_y_to_x(this%rbuffy1,this%rbuffx1)
    call this%get_ddx(this%rbuffx1,this%rbuffx2)
    call transpose_x_to_y(this%rbuffx2,this%rbuffy2,this%gp)
    call transpose_y_to_z(this%rbuffy2,dudx,this%gp)

end subroutine 


subroutine get_ddz(this,u,dudz)
    class(d3q19), intent(in) :: this
    real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)), intent(in ) :: u
    real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)), intent(out) :: dudz
    integer :: nz

    nz = this%gp%zsz(3)
    
    ! dudz
    if (this%iszPeriodic) then
        dudz(:,:,2:nz-1) = half*(u(:,:,3:nz) - u(:,:,1:nz-2))
        dudz(:,:,1)      = half*(u(:,:,2)    - u(:,:,nz)    )
        dudz(:,:,nz)     = half*(u(:,:,1)    - u(:,:,nz-1)  )
    else
        dudz(:,:,2:nz-1) = half*(u(:,:,3:nz) - u(:,:,1:nz-2))
        !dudz(:,:,1) = (-3.d0/2.d0)*u(:,:,1) + 2.d0*u(:,:,2) - 0.5d0*(u(:,:,3)) 
        !dudz(:,:,nz) = -((-3.d0/2.d0)*u(:,:,nz) + 2.d0*u(:,:,nz-1) - 0.5d0*(u(:,:,nz-2)))
        dudz(:,:,1)  =   (-11.d0/6.d0)*u(:,:,1 ) + 3.d0*u(:,:,2   ) - (3.d0/2.d0)*u(:,:,3   ) + (1.d0/3.d0)*u(:,:,4   )
        dudz(:,:,nz) = -((-11.d0/6.d0)*u(:,:,nz) + 3.d0*u(:,:,nz-1) - (3.d0/2.d0)*u(:,:,nz-2) + (1.d0/3.d0)*u(:,:,nz-3))
    end if 

end subroutine


subroutine get_ddy(this,u,dudy)
    class(d3q19), intent(in) :: this
    real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3)), intent(in ) :: u
    real(rkind), dimension(this%gp%ysz(1),this%gp%ysz(2),this%gp%ysz(3)), intent(out) :: dudy
    integer :: ny 

    ny = this%gp%ysz(2)

    dudy(:,2:ny-1,:) = half*(u(:,3:ny,:) - u(:,1:ny-2,:))
    dudy(:,1,:) = half*(u(:,2,:) - u(:,ny,:))
    dudy(:,ny,:) = half*(u(:,1,:) - u(:,ny-1,:))
    
end subroutine

subroutine get_ddx(this,u,dudx)
    class(d3q19), intent(in) :: this
    real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(in ) :: u
    real(rkind), dimension(this%gp%xsz(1),this%gp%xsz(2),this%gp%xsz(3)), intent(out) :: dudx
    integer :: nx 

    nx = this%gp%xsz(1)

    dudx(2:nx-1,:,:) = half*(u(3:nx,:,:) - u(1:nx-2,:,:))
    dudx(1,:,:) = half*(u(2,:,:) - u(nx,:,:))
    dudx(nx,:,:) = half*(u(1,:,:) - u(nx-1,:,:))
end subroutine
