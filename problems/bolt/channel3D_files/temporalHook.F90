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
    real(rkind) :: fmin, fmax
    real(rkind) :: utau_mn, ubc_mn

    utau_mn = 0.5d0*(p_sum(utau_up) + p_sum(utau_do))/real(bgp%nx*bgp%ny,rkind)
    ubc_mn = 0.5d0*(p_sum(ubc_up) + p_sum(ubc_do))/real(bgp%nx*bgp%ny,rkind)

    fmax = p_maxval(maxval(bgp%f))
    fmin = p_minval(minval(bgp%f))
    call message(0,"Time:", bgp%getPhysTime())
    call message(0,"TIDX:", bgp%step)
    call message(1,"u_tau:", utau_mn)
    call message(1,"u_bc:", ubc_mn)
    call message_min_max(1,"Bounds for f  :", fmin, fmax) 
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
