module temporalHook
    use kind_parameters, only: rkind
    use d3q19mod, only: d3q19
    use exits, only: message, message_min_max
    use timer, only: tic, toc
    use reductions 
    implicit none
contains 

subroutine doTemporalStuff(bgp)
    class(d3q19), intent(in) :: bgp
    real(rkind) :: fmin, fmax


    fmax = p_maxval(maxval(bgp%f))
    fmin = p_minval(minval(bgp%f))
    
    
    call message(0,"Time:", bgp%getPhysTime())
    call message(0,"TIDX:", bgp%step)
    call message_min_max(1,"Bounds for f  :", fmin, fmax) 
    call message_min_max(1,"Bounds for u  :", p_minval(minval(bgp%ux))*bgp%delta_u, p_maxval(maxval(bgp%ux))*bgp%delta_u)
    call message_min_max(1,"Bounds for v  :", p_minval(minval(bgp%uy))*bgp%delta_u, p_maxval(maxval(bgp%uy))*bgp%delta_u)
    call message_min_max(1,"Bounds for w  :", p_minval(minval(bgp%uz))*bgp%delta_u, p_maxval(maxval(bgp%uz))*bgp%delta_u)
    call message_min_max(1,"Bounds for rho:", p_minval(minval(bgp%rho)), p_maxval(maxval(bgp%rho)))
    call toc()
    call message("========================================")
    call tic()
end subroutine 
end module 
