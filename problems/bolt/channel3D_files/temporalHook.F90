module temporalHook
    use kind_parameters, only: rkind
    use d3q19mod, only: d3q19
    use exits, only: message, message_min_max
    use d3q19_channel3D
    use timer, only: tic, toc
    use reductions 
    implicit none
contains 

subroutine doTemporalStuff(bgp)
    class(d3q19), intent(in) :: bgp
    real(rkind) :: usum, vsum, wsum

    usum = p_sum(sum(bgp%ux))/real(bgp%nx*bgp%ny*bgp%nz,rkind) 
    vsum = p_sum(sum(bgp%uy))/real(bgp%nx*bgp%ny*bgp%nz,rkind)
    wsum = p_sum(sum(bgp%uz))/real(bgp%nx*bgp%ny*bgp%nz,rkind)

    call message(0,"Time:", bgp%getPhysTime())
    call message(0,"TIDX:", bgp%step)
    call message(1,"u_tau:", utau)
    call message(1,"u_match:", umatch)
    call message(1,"Domain u_sum:", usum)
    call message(1,"Domain v_sum:", vsum)
    call message(1,"Domain w_sum:", wsum)
    call message_min_max(1,"Bounds for u  :", p_minval(minval(bgp%ux))*bgp%delta_u, p_maxval(maxval(bgp%ux))*bgp%delta_u)
    call message_min_max(1,"Bounds for v  :", p_minval(minval(bgp%uy))*bgp%delta_u, p_maxval(maxval(bgp%uy))*bgp%delta_u)
    call message_min_max(1,"Bounds for w  :", p_minval(minval(bgp%uz))*bgp%delta_u, p_maxval(maxval(bgp%uz))*bgp%delta_u)
    call message_min_max(1,"Bounds for rho:", p_minval(minval(bgp%rho)), p_maxval(maxval(bgp%rho)))
    call message_min_max(1,"Bounds for nuT:", p_minval(minval(bgp%nuSGS)), p_maxval(maxval(bgp%nuSGS)))
    call toc()
    call message("========================================")
    call tic()
end subroutine 
end module 
