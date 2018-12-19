module temporalHook
    use kind_parameters, only: rkind
    use d3q19mod, only: d3q19
    use exits, only: message, message_min_max
    use d3q19_taylorgreen3D
    use timer, only: tic, toc
    use reductions 
    implicit none
contains 

subroutine doTemporalStuff(bgp)
    class(d3q19), intent(inout) :: bgp
    real(rkind) :: fmin, fmax, DomainTKE
    real(rkind) :: error_u, error_v, error_w

    fmax = p_maxval(maxval(bgp%f))
    fmin = p_minval(minval(bgp%f))
    call message(0,"Time:", bgp%getPhysTime())
    call message(0,"TIDX:", bgp%step)
    call message_min_max(1,"Bounds for f  :", fmin, fmax) 
    call message_min_max(1,"Bounds for u  :", p_minval(minval(bgp%ux))*bgp%delta_u, p_maxval(maxval(bgp%ux))*bgp%delta_u)
    call message_min_max(1,"Bounds for v  :", p_minval(minval(bgp%uy))*bgp%delta_u, p_maxval(maxval(bgp%uy))*bgp%delta_u)
    call message_min_max(1,"Bounds for w  :", p_minval(minval(bgp%uz))*bgp%delta_u, p_maxval(maxval(bgp%uz))*bgp%delta_u)
    call message_min_max(1,"Bounds for rho:", p_minval(minval(bgp%rho)), p_maxval(maxval(bgp%rho)))
    call message("----------------------------------------")
    bgp%buff1 = bgp%ux*bgp%ux + bgp%uy*bgp%uy + bgp%uz*bgp%uz
    DomainTKE = 0.5d0*bgp%delta_u*bgp%delta_u*p_sum(sum(bgp%buff1))/(bgp%nx*bgp%ny*bgp%nz)
    call message(1,"Domain TKE:", DomainTKE)
    call message("----------------------------------------")
    call toc()
    call message("========================================")
    call tic()
end subroutine 
end module 
