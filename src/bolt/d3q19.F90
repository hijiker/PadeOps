module d3q19mod
    use kind_parameters, only: rkind, clen 
    use exits, only: GracefulExit, message, check_exit
    use decomp_2d

    implicit none

    private

    public :: d3q19 


    type :: d3q19
        private
        
        type(decomp_info) :: gp
        integer :: nx, ny, nz

        character(len=clen) :: inputdir, outputdir

        real(rkind), dimension(:,:,:,:), allocatable :: fnew, fstar, feq, fneq

        real(rkind) :: delta_t, delta_x, delta_u
        real(rkind) :: Re

        real(rkind), dimension(:,:,:), allocatable :: u, v, w, nu, tau

        contains
            procedure :: init
            procedure :: destroy

            procedure :: collide 
            procedure :: stream

    end type


contains

    subroutine init(this,inputfile)
        class(d3q19), intent(inout) :: this
        character(len=*), intent(in) :: inputfile 




    end subroutine 

    subroutine destroy(this)
        class(d3q19), intent(inout) :: this


    end subroutine 


    subroutine collide(this)
        class(d3q19), intent(inout) :: this


    end subroutine


    subroutine stream(this)
        class(d3q19), intent(inout) :: this

    end subroutine 

end module 
