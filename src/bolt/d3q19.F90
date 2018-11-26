module d3q19mod
    use kind_parameters, only: rkind, clen 
    use exits, only: GracefulExit, message, check_exit
    use decomp_2d
    use constants, only: one, zero, two, half 
    use bolt_hooks
    implicit none

    private
    public :: d3q19, test_lattice_definition 

    integer, parameter :: nvels = 19
    real(rkind), parameter :: o18 = 1._rkind/18._rkind
    real(rkind), parameter :: o3  = 1._rkind/3._rkind
    real(rkind), parameter :: o36 = 1._rkind/36._rkind

    real(rkind), parameter :: pon  =  one
    real(rkind), parameter :: mon  = -one
    real(rkind), parameter :: zer  =  zero 


    real(rkind), dimension(nvels), parameter :: w  = [o18,o18,o18,o18,o18,o18,o36,o36,o36,o36,o36,o36,o36,o36,o36,o36,o36,o36,o3 ]
    real(rkind), dimension(nvels), parameter :: cx = [pon,mon,zer,zer,zer,zer,pon,mon,pon,mon,zer,zer,pon,mon,pon,mon,zer,zer,zer]
    real(rkind), dimension(nvels), parameter :: cy = [zer,zer,pon,mon,zer,zer,pon,mon,zer,zer,pon,mon,mon,pon,zer,zer,pon,mon,zer]
    real(rkind), dimension(nvels), parameter :: cz = [zer,zer,zer,zer,pon,mon,zer,zer,pon,mon,pon,mon,zer,zer,mon,pon,mon,pon,zer]

    real(rkind), parameter :: csq = 1._rkind/3._rkind
    real(rkind), parameter :: onebycsq = 3._rkind
    real(rkind), parameter :: oneby2c4 = 9._rkind/2._rkind

    integer, dimension(9), parameter :: fparallel = [1,2,3,4,7,8,13,14,19]
    integer, dimension(5), parameter :: fplus     = [5,9,11,16,18]
    integer, dimension(5), parameter :: fminus    = [6,10,12,15,17]

    type :: d3q19

        type(decomp_info), public :: gp
        integer :: nx, ny, nz

        character(len=clen) :: inputdir, outputdir

        real(rkind), dimension(:,:,:,:), allocatable :: f

        real(rkind) :: delta_t, delta_x, delta_u, delta_nu
        real(rkind) :: Re

        real(rkind), dimension(:,:,:), allocatable :: rho, ux, uy, uz, nu, tau, VisBuff
        real(rkind) :: Fx = zero, Fy = zero, Fz = zero 
        real(rkind), dimension(:,:,:), allocatable :: OneByTau, OneByTwoTau
        real(rkind), dimension(3,3,nvels) :: QTensor
        real(rkind), dimension(:,:,:), allocatable :: PiTensor

        logical :: useRestart = .false., useConstantBodyForce, isZPeriodic = .true., useSmagorinsky=.false.

        integer :: CollisionModel, restart_RunID, restart_timeID, RunID, tid_vis=10000, tid_restart=10000

        integer :: step 

        ! Boundary condition data
        real(rkind), dimension(:,:), allocatable :: rhoBC, uTop, vTop, wTop, uBot, vBot, wBot
        real(rkind), dimension(:,:,:), allocatable :: feqBC, fneqBC, PiBC

        ! Streaming data arrays
        real(rkind), dimension(:,:), allocatable :: YZslice, XZslice, XYslice
        real(rkind), dimension(:,:), allocatable :: YZsendb, XZsendb, XYsendb
        real(rkind), dimension(:), allocatable :: z1d_send, z1d_recv
        integer :: XneighLeft, XneighRight, YneighDown, YneighUp

        contains
            procedure :: init
            procedure :: destroy
            procedure :: time_advance
            procedure :: collide 
            procedure :: stream
            procedure :: compute_macroscopic 
            procedure :: updateBCs
            procedure :: GetPhysTime

            procedure :: dumpVisualizationFields
            procedure :: wrapup_timestep

            procedure, private :: collision_BGK 
            procedure, private :: collision_BGK_Force 
            procedure, private :: NEQ_BB_Reg
            procedure, private :: compute_rhoBC
            procedure, private :: RegularizeFneq_BC
            procedure, private :: compute_f_BC
            procedure, private :: initialize_f_to_feq 
            procedure, private :: dumpRestart
            procedure, private :: readRestart
            procedure, private :: dumpVisBuff
            procedure, private :: get_processor_topo
            procedure, private :: allocate_lattice_memory
            procedure, private :: read_inputfile
    end type


contains

#include "domaindecomp_code/Decompose_LBM_z.F90"
#include "domaindecomp_code/allocate_deallocate.F90"

#include "equilibrium_code/second_order_eq.F90"
#include "testing_code/test_lattice_def.F90"
#include "io_code/all_io.F90"

#include "collision_code/BGK.F90"

#include "conversions_code/f_to_macro.F90"
#include "conversions_code/macro_to_feq.F90"

#include "d3q19_codes/d3q19_streaming.F90"
#include "d3q19_codes/d3q19_init.F90"

#include "BC_code/NEQ_BB_reg.F90"

    subroutine time_advance(this)
        class(d3q19), intent(inout) :: this

        call this%collide()
        call this%stream()   ! Look at "d3q19_codes/d3q19_streaming.F90"
        call this%updateBCs()
        call this%wrapup_timestep()

    end subroutine 

    pure function getPhysTime(this) result(time)
        class(d3q19), intent(in) :: this
        real(rkind) :: time

        time = this%step*this%delta_t
    end function 

    subroutine updateBCs(this)
        class(d3q19), intent(inout) :: this
        
        if (this%isZperiodic) then
            return 
        else
            ! User defnined BC (dirichlet)
            call getWallBC_bolt(this%gp, this%ux, this%uy, this%uz, this%uBot, &
                & this%vBot, this%wBot, this%uTop, this%vTop, this%wTop)

            call this%NEQ_BB_Reg()
        end if 

    end subroutine 

    subroutine wrapup_timestep(this)
        class(d3q19), intent(inout) :: this

        this%step = this%step + 1
        call this%compute_macroscopic()
        
        if (mod(this%step,this%tid_vis) == 0) then
            call this%dumpVisualizationFields()
        end if 

        if (mod(this%step,this%tid_restart)==0) then
            call this%dumpRestart()
        end if 

    end subroutine 

    subroutine collide(this)
        class(d3q19), intent(inout) :: this

        select case(this%CollisionModel) 
        case (0) ! BGK
            if (this%useConstantBodyForce) then
                call this%collision_BGK_Force()
            else
                call this%collision_BGK()
            end if 
        case default
            call this%collision_BGK_Force()
        end select 

    end subroutine 
   
end module 
