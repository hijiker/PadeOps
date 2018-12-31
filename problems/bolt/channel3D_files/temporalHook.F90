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
    real(rkind) :: utau_mn

    fmax = p_maxval(maxval(bgp%f))
    fmin = p_minval(minval(bgp%f))
    
    
    call message(0,"Time:", bgp%getPhysTime())
    call message(0,"TIDX:", bgp%step)
    call message(0,"u_mean:", umean)
    call message_min_max(1,"Bounds for f  :", fmin, fmax) 
    call message_min_max(1,"Bounds for u  :", p_minval(minval(bgp%ux))*bgp%delta_u, p_maxval(maxval(bgp%ux))*bgp%delta_u)
    call message_min_max(1,"Bounds for v  :", p_minval(minval(bgp%uy))*bgp%delta_u, p_maxval(maxval(bgp%uy))*bgp%delta_u)
    call message_min_max(1,"Bounds for w  :", p_minval(minval(bgp%uz))*bgp%delta_u, p_maxval(maxval(bgp%uz))*bgp%delta_u)
    call message_min_max(1,"Bounds for rho:", p_minval(minval(bgp%rho)), p_maxval(maxval(bgp%rho)))
    if (bgp%useSGSmodel) then
        call message_min_max(1,"Bounds for nut:", bgp%delta_nu*p_minval(minval(bgp%nuSGS)), bgp%delta_nu*p_maxval(maxval(bgp%nuSGS)))
    end if 
    if ((allocated(utau_up)) .and. (allocated(utau_do))) then
        utau_mn = 0.5*(p_sum(sum(utau_up)) + p_sum(sum(utau_do)))/real(bgp%gp%xsz(1)*bgp%gp%ysz(2),rkind)
        call message(1, "u_tau", utau_mn)
    end if 
    call toc()
    call message("========================================")
    call tic()
end subroutine 
end module 
