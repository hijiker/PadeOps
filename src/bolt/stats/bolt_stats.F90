subroutine start_stats(this)
    class(d3q19), intent(inout) :: this
    integer, parameter :: n_means = 11

    allocate(this%dat_array(this%nz,n_means))
    allocate(this%sum_array(this%nz,n_means))
    allocate(this%stats_1d(this%nz))
    this%dat_count = 0
    this%sum_array = 0.d0 

end subroutine 

subroutine end_stats(this)
    class(d3q19), intent(inout) :: this

    deallocate(this%dat_array)
    deallocate(this%sum_array)
    deallocate(this%stats_1d)
end subroutine


subroutine calculate_stats(this)
    class(d3q19), intent(inout) :: this

    call this%average_xy(this%ux)
    this%sum_array(:,1) = this%sum_array(:,1) + this%stats_1d
    
    call this%average_xy(this%uy)
    this%sum_array(:,2) = this%sum_array(:,2) + this%stats_1d
    
    call this%average_xy(this%uz)
    this%sum_array(:,3) = this%sum_array(:,3) + this%stats_1d
   
    this%rbuffz1 = this%ux*this%ux
    call this%average_xy(this%rbuffz1)
    this%sum_array(:,4) = this%sum_array(:,4) + this%stats_1d
    
    this%rbuffz1 = this%ux*this%uy
    call this%average_xy(this%rbuffz1)
    this%sum_array(:,5) = this%sum_array(:,5) + this%stats_1d
    
    this%rbuffz1 = this%ux*this%uz
    call this%average_xy(this%rbuffz1)
    this%sum_array(:,6) = this%sum_array(:,6) + this%stats_1d
    
    this%rbuffz1 = this%uy*this%uy
    call this%average_xy(this%rbuffz1)
    this%sum_array(:,7) = this%sum_array(:,7) + this%stats_1d
    
    this%rbuffz1 = this%uy*this%uz
    call this%average_xy(this%rbuffz1)
    this%sum_array(:,8) = this%sum_array(:,8) + this%stats_1d

    this%rbuffz1 = this%uz*this%uz
    call this%average_xy(this%rbuffz1)
    this%sum_array(:,9) = this%sum_array(:,9) + this%stats_1d
    
    call this%average_xy(this%rho)
    this%sum_array(:,10) = this%sum_array(:,10) + this%stats_1d

    call this%average_xy(this%nusgs)
    this%sum_array(:,11) = this%sum_array(:,11) + this%stats_1d

    this%dat_count = this%dat_count + 1

end subroutine 


subroutine average_xy(this,f,ffluct)
    use reductions, only: p_sum
    class(d3q19), intent(inout) :: this
    real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)), intent(in) :: f
    real(rkind), dimension(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)), intent(out), optional :: ffluct
    integer :: k 
    real(rkind) :: mfact

    mfact = 1.d0/real(this%nx*this%ny,rkind)

    do k = 1,this%nz
        this%stats_1d(k) = p_sum(sum(f(:,:,k)))*mfact
    end do 

    if (present(ffluct)) then
        do k = 1,this%nz
            ffluct(:,:,k) = f(:,:,k) - this%stats_1d(k)
        end do 
    end if 

end subroutine 
