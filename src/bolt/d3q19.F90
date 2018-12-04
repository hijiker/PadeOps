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

    real(rkind), parameter :: C_smag = 0.16_rkind 
    real(rkind), parameter :: C_smag_sq = C_smag*C_smag

    type :: d3q19

        type(decomp_info), public :: gp
        integer :: nx, ny, nz

        character(len=clen) :: inputdir, outputdir

        real(rkind), dimension(:,:,:,:), allocatable :: f

        real(rkind) :: delta_t, delta_x, delta_u, delta_nu
        real(rkind) :: Re

        real(rkind), dimension(:,:,:), allocatable :: rho, ux, uy, uz,  tau, VisBuff
        real(rkind), dimension(:,:,:), allocatable :: Fx, Fy, Fz
        real(rkind) :: nu
        real(rkind), dimension(:,:,:), allocatable :: OneByTau, OneByTwoTau
        real(rkind), dimension(3,3,nvels) :: QTensor

        logical :: useRestart = .false., useConstantBodyForce = .false., isZPeriodic = .true.
        logical :: restartWithTau = .false. 

        integer :: CollisionModel, restart_RunID, restart_timeID, RunID, tid_vis=10000, tid_restart=10000

        integer :: step 

        logical :: useSGSmodel = .false., useSpaceTimeBodyForce = .false. 
 
        ! Boundary condition data
        real(rkind), dimension(:,:), allocatable :: rhoBC, uTop, vTop, wTop, uBot, vBot, wBot
        real(rkind), dimension(:,:), allocatable :: tau_T, tau_B
        real(rkind), dimension(:,:,:), allocatable :: feqBC, fneqBC
        real(rkind), dimension(:,:,:,:), allocatable :: PiBC
        real(rkind), dimension(:,:,:,:), allocatable :: PiTensor
        real(rkind), dimension(:,:,:), allocatable :: buff1, buff2, nuSGS

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

            procedure :: compute_Pi
            procedure :: dumpVisualizationFields
            procedure :: wrapup_timestep
            procedure :: update_bodyForce

            procedure, private :: compute_tau_smag
            procedure, private :: collision_BGK 
            procedure, private :: collision_BGK_Force 
            procedure, private :: RegBGK_Force 
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

! GENERIC PROCEDURES (FOR AN ARBITRARY LATTICE)
#include "domaindecomp_code/Decompose_LBM_z.F90"
#include "domaindecomp_code/allocate_deallocate.F90"
#include "equilibrium_code/second_order_eq.F90"
#include "testing_code/test_lattice_def.F90"
#include "io_code/all_io.F90"
#include "collision_code/BGK.F90"
#include "collision_code/RegBGK.F90"
#include "conversions_code/f_to_macro.F90"
#include "conversions_code/macro_to_feq.F90"
#include "BC_code/NEQ_BB_reg.F90"
#include "SGS_code/Smag.F90"

! D3Q19 SPECIFIC PROCEDURES
#include "d3q19_codes/d3q19_streaming.F90"
#include "d3q19_codes/d3q19_init.F90"

    subroutine time_advance(this)
        use timer, only: tic, toc
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
            call getWallBC_bolt(this%gp, this%Re, this%delta_u, this%ux, this%uy, this%uz, this%uBot, &
                & this%vBot, this%wBot, this%uTop, this%vTop, this%wTop)

            call this%NEQ_BB_Reg()
        end if 

    end subroutine 

    subroutine wrapup_timestep(this)
        class(d3q19), intent(inout) :: this

        this%step = this%step + 1

        if (mod(this%step,this%tid_restart)==0) then
            call this%dumpRestart()
        end if 

        if (this%useSGSmodel) then
            call this%compute_tau_smag()
        end if

        if (this%useSpaceTimeBodyForce) then
            call getBodyForce(this%gp, this%getPhysTime(),this%delta_u, this%delta_t, &
              &  this%ux, this%uy, this%uz, this%Fx, this%Fy, this%Fz)
        end if 
        
        if (mod(this%step,this%tid_vis) == 0) then
            call this%dumpVisualizationFields()
        end if 
        
        call this%compute_macroscopic()
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
        case (1) ! Reg. BGK
            call this%RegBGK_Force()
        case default
            call this%collision_BGK_Force()
        end select 

    end subroutine 
   
    subroutine update_bodyForce(this,ForceX)
        class(d3q19), intent(inout) :: this
        real(rkind), intent(in) :: ForceX

        this%Fx = ForceX*this%delta_t/this%delta_u

    end subroutine 

    subroutine compute_Pi(this)
        class(d3q19), intent(inout) :: this
        real(rkind), dimension(nvels) :: feq, Fvals
        integer :: i, j, k, idx
        real(rkind) :: fneq

        do k = 1,this%gp%zsz(3)
            do j = 1,this%gp%zsz(2)
                do i = 1,this%gp%zsz(1)
                    this%Pitensor(:,i,j,k) = zero
                    do idx = 1,nvels
                        call get_Feq_2ndOrder(this%ux(i,j,k),this%uy(i,j,k),this%uz(i,j,k), &
                                & this%rho(i,j,k),idx,this%Qtensor(:,:,idx),feq(idx))
                        
                        call get_ForceSource_2ndOrder(this%ux(i,j,k),this%uy(i,j,k),this%uz(i,j,k), &
                                & this%Fx(i,j,k), this%Fy(i,j,k), this%Fz(i,j,k), idx, &
                                & this%Qtensor(:,:,idx), Fvals(idx))
                        
                        fneq = this%f(i,j,k,idx) - feq(idx)
                        
                        this%Pitensor(1,i,j,k) = this%Pitensor(1,i,j,k) + cx(idx)*cx(idx)*(fneq + half*Fvals(idx))
                        this%Pitensor(2,i,j,k) = this%Pitensor(2,i,j,k) + cx(idx)*cy(idx)*(fneq + half*Fvals(idx))
                        this%Pitensor(3,i,j,k) = this%Pitensor(3,i,j,k) + cx(idx)*cz(idx)*(fneq + half*Fvals(idx))
                        this%Pitensor(4,i,j,k) = this%Pitensor(4,i,j,k) + cy(idx)*cy(idx)*(fneq + half*Fvals(idx))
                        this%Pitensor(5,i,j,k) = this%Pitensor(5,i,j,k) + cy(idx)*cz(idx)*(fneq + half*Fvals(idx))
                        this%Pitensor(6,i,j,k) = this%Pitensor(6,i,j,k) + cz(idx)*cz(idx)*(fneq + half*Fvals(idx))
                    end do
                end do
            end do 
        end do 
    
    end subroutine 


end module 
