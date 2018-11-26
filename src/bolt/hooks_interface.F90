module bolt_hooks
    use kind_parameters, only: rkind, clen
    use decomp_2d, only: decomp_info
    implicit none

    interface initfields_bolt
        subroutine initfields_bolt(decomp, inputfile, delta_x, rho, ux, uy, uz)
            import :: rkind
            import :: clen
            import :: decomp_info
            type(decomp_info), intent(in) :: decomp
            character(len=clen), intent(in) :: inputfile
            real(rkind), intent(in) :: delta_x
            real(rkind), dimension(:,:,:), intent(out) :: rho, ux, uy, uz  

        end subroutine 
    end interface

    interface getWallBC_bolt
        subroutine getWallBC_bolt(decomp, ux, uy, uz, uxB, uyB, uzB, uxT, uyT, uzT)
            import :: rkind
            import :: decomp_info
            type(decomp_info), intent(in) :: decomp
            real(rkind), dimension(:,:,:), intent(in) :: ux, uy, uz
            real(rkind), dimension(:,:), intent(out) :: uxB,uyB,uzB,uxT,uyT,uzT
        end subroutine 
    end interface

end module 
