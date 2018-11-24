module d3q19mod
    use kind_parameters, only: rkind, clen 
    use exits, only: GracefulExit, message, check_exit
    use decomp_2d
    use constants, only: one, zero, two, half 
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

    integer, dimension(nvels), parameter :: oppo = [2,1,4,3,6,5,8,9,10,9,12,11,14,13,16,15,18,17,19]

    type :: d3q19
        private

        type(decomp_info), public :: gp
        integer, public :: nx, ny, nz

        character(len=clen) :: inputdir, outputdir

        real(rkind), dimension(:,:,:,:), allocatable, public :: f

        real(rkind) :: delta_t, delta_x, delta_u
        real(rkind) :: Re

        real(rkind), dimension(:,:,:), allocatable :: rho, ux, uy, uz, nu, tau, Fx, Fy, Fz
        real(rkind), dimension(:,:,:), allocatable :: OneByTau, OneByTwoTau
        real(rkind), dimension(3,3,nvels) :: QTensor
        real(rkind), dimension(:,:,:), allocatable :: PiTensor

        logical :: useBodyForce 

        real(rkind), dimension(:,:), allocatable :: YZslice, XZslice, XYslice
        real(rkind), dimension(:,:), allocatable :: YZsendb, XZsendb, XYsendb
        real(rkind), dimension(:), allocatable :: z1d_send, z1d_recv

        integer :: XneighLeft, XneighRight, YneighDown, YneighUp
        contains
            procedure :: init
            procedure :: destroy
            procedure :: collide 
            procedure :: stream
            procedure :: compute_macroscopic 

            procedure, private :: get_processor_topo
            procedure, private :: allocate_lattice_memory

    end type


contains

#include "streaming_code/d3q19_streaming.F90"
#include "domaindecomp_code/Decompose_LBM_z.F90"


    subroutine init(this,inputfile)
        class(d3q19), intent(inout) :: this
        character(len=*), intent(in) :: inputfile 
        integer :: idx

        ! Populate the Qtensor
        do idx = 1,nvels
            this%Qtensor(1,1,idx) = cx(idx)*cx(idx) - csq
            this%Qtensor(2,2,idx) = cy(idx)*cy(idx) - csq
            this%Qtensor(3,3,idx) = cz(idx)*cz(idx) - csq

            this%Qtensor(1,2,idx) = cx(idx)*cy(idx) 
            this%Qtensor(2,1,idx) = cx(idx)*cy(idx) 

            this%Qtensor(1,3,idx) = cx(idx)*cz(idx) 
            this%Qtensor(3,1,idx) = cx(idx)*cz(idx) 

            this%Qtensor(2,3,idx) = cy(idx)*cz(idx) 
            this%Qtensor(3,2,idx) = cy(idx)*cz(idx) 
        end do 


        call this%get_processor_topo()
        call this%allocate_lattice_memory()

    end subroutine 

    subroutine allocate_lattice_memory(this)
        class(d3q19), intent(inout) :: this

        allocate(this%f  (this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3),nvels))
        allocate(this%ux (this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        allocate(this%uy (this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        allocate(this%uz (this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        allocate(this%rho(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        
        allocate(this%Fx(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        allocate(this%Fy(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))
        allocate(this%Fz(this%gp%zsz(1),this%gp%zsz(2),this%gp%zsz(3)))

        allocate(this%XYslice(this%gp%zsz(1),this%gp%zsz(2)))
        allocate(this%XZslice(this%gp%zsz(1),this%gp%zsz(3)))
        allocate(this%YZslice(this%gp%zsz(2),this%gp%zsz(3)))
        
        allocate(this%YZsendb(this%gp%zsz(2),this%gp%zsz(3)))
        allocate(this%XZsendb(this%gp%zsz(1),this%gp%zsz(3)))
        allocate(this%z1d_send(this%gp%zsz(3)))
        allocate(this%z1d_recv(this%gp%zsz(3)))

    end subroutine 
    
    subroutine destroy(this)
        class(d3q19), intent(inout) :: this

        deallocate(this%f, this%ux, this%uy, this%uz, this%rho)
        deallocate(this%XYslice, this%XZslice, this%YZslice)
        deallocate(this%YZsendb, this%XZsendb)
        deallocate(this%z1d_send, this%z1d_recv)

    end subroutine 


    subroutine collide(this)
        class(d3q19), intent(inout) :: this
        integer :: i, j, k, idx 
        real(rkind) :: OneByTau, feq, Force

        do idx = 1,nvels
            do k = 1,this%gp%zsz(3)
                do j = 1,this%gp%zsz(2)
                    !$omp simd 
                    do i = 1,this%gp%zsz(1)
                        call get_Feq_2ndOrder(this%ux(i,j,k),this%uy(i,j,k),this%uz(i,j,k), &
                                & this%rho(i,j,k),idx,this%Qtensor(:,:,idx),feq)
                        
                        call get_ForceSource_2ndOrder(this%ux(i,j,k),this%uy(i,j,k),this%uz(i,j,k), &
                                & this%Fx(i,j,k), this%Fy(i,j,k), this%Fz(i,j,k), idx, &
                                & this%Qtensor(:,:,idx), Force)
                        
                        oneByTau = one/this%tau(i,j,k)
                        
                        this%f(i,j,k,idx) = (one - oneByTau)*this%f(i,j,k,idx) + oneByTau*feq &
                                        & + (one - half*oneBytau)*Force
                    end do 
                end do 
            end do 
        end do 

    end subroutine

    pure subroutine get_Feq_2ndOrder(ux,uy,uz,rho,fidx,Qtensor,feq)
        real(rkind), intent(in) :: ux, uy, uz, rho
        integer, intent(in) :: fidx
        real(rkind), dimension(3,3), intent(in) :: Qtensor
        real(rkind), intent(out) :: feq

        real(rkind) :: first, second

        first = cx(fidx)*ux + cy(fidx)*uy + cz(fidx)*uz

        second =     ux*ux*Qtensor(1,1) +     uy*uy*Qtensor(2,2) +     uz*uz*Qtensor(3,3) &
             & + two*ux*uy*Qtensor(1,2) + two*ux*uz*Qtensor(1,3) + two*uy*uz*Qtensor(2,3) 
        
        feq = w(fidx)*rho*(one + onebycsq*first + oneby2c4*second)

    end subroutine 

    pure subroutine get_ForceSource_2ndOrder(ux,uy,uz,ForceX,ForceY,ForceZ,fidx,Qtensor,Fsource)
        real(rkind), intent(in) :: ux, uy, uz, ForceX, ForceY, ForceZ
        integer, intent(in) :: fidx
        real(rkind), dimension(3,3), intent(in) :: Qtensor
        real(rkind), intent(out) :: Fsource

        real(rkind) :: first, second

        first = cx(fidx)*ForceX + cy(fidx)*ForceY + cz(fidx)*ForceZ

        second = (Qtensor(1,1)*ux + Qtensor(1,2)*uy + Qtensor(1,3)*uz)*ForceX &
               + (Qtensor(2,1)*ux + Qtensor(2,2)*uy + Qtensor(2,3)*uz)*ForceY &
               + (Qtensor(3,1)*ux + Qtensor(3,2)*uy + Qtensor(3,3)*uz)*ForceZ
        
        Fsource = w(fidx)*(onebycsq*first + oneby2c4*second)

    end subroutine 



    subroutine compute_macroscopic(this)
        class(d3q19), intent(inout) :: this
        integer :: i, j, k, onebyrho

        do k = 1,this%gp%zsz(3)
            do j = 1,this%gp%zsz(2)
                !$omp simd
                do i = 1,this%gp%zsz(1)
                    this%rho(i,j,k) = this%f(i,j,k,1 ) + this%f(i,j,k,2 ) + this%f(i,j,k,3 ) + this%f(i,j,k,4 ) &
                                 &  + this%f(i,j,k,5 ) + this%f(i,j,k,6 ) + this%f(i,j,k,7 ) + this%f(i,j,k,8 ) &
                                 &  + this%f(i,j,k,9 ) + this%f(i,j,k,10) + this%f(i,j,k,11) + this%f(i,j,k,12) &
                                 &  + this%f(i,j,k,13) + this%f(i,j,k,14) + this%f(i,j,k,15) + this%f(i,j,k,16) &
                                 &  + this%f(i,j,k,17) + this%f(i,j,k,18) + this%f(i,j,k,19)

                    onebyrho = one/this%rho(i,j,k)  

                    this%ux(i,j,k)  = (this%f(i,j,k,1 ) - this%f(i,j,k,2 ) + this%f(i,j,k,7 ) - this%f(i,j,k,8 ) &
                                 &  +  this%f(i,j,k,9 ) - this%f(i,j,k,10) + this%f(i,j,k,13) - this%f(i,j,k,14) &
                                 &  +  this%f(i,j,k,15) - this%f(i,j,k,16))*onebyrho                           
                                                                                                               
                    this%uy(i,j,k)  = (this%f(i,j,k,3 ) - this%f(i,j,k,4 ) + this%f(i,j,k,7 ) - this%f(i,j,k,8 ) &
                                 &  +  this%f(i,j,k,11) - this%f(i,j,k,12) - this%f(i,j,k,13) + this%f(i,j,k,14) &
                                 &  +  this%f(i,j,k,17) - this%f(i,j,k,18))*onebyrho                           
                                                                                                               
                    this%uz(i,j,k)  = (this%f(i,j,k,5 ) - this%f(i,j,k,6 ) + this%f(i,j,k,9 ) - this%f(i,j,k,10) &
                                 &  +  this%f(i,j,k,11) - this%f(i,j,k,12) - this%f(i,j,k,15) + this%f(i,j,k,16) &
                                 &  -  this%f(i,j,k,17) + this%f(i,j,k,18))*onebyrho 
               
                end do 
            end do 
        end do 

    end subroutine 


    subroutine test_lattice_definition()
        real(rkind) :: s

        call message(0,"PERFORMING LATTICE TESTS")
        ! Test 1: sum of weights
        s = sum(w); call message(1,"TEST 1: Sum of weights ",s)

        ! Test 2: wi ci
        s = sum(w*cx); call message(1,"TEST 2a: Sum of {w_i c_i_x} ",s)
        s = sum(w*cy); call message(1,"TEST 2b: Sum of {w_i c_i_y} ",s)
        s = sum(w*cz); call message(1,"TEST 2c: Sum of {w_i c_i_z} ",s)

        ! Test 3: wi ci_alpha ci_beta
        s = sum(w*cx*cx); call message(1,"TEST 3a: Sum of {w_i c_i_x c_i_x} ",s)
        s = sum(w*cx*cy); call message(1,"TEST 3a: Sum of {w_i c_i_x c_i_y} ",s)
        s = sum(w*cx*cz); call message(1,"TEST 3a: Sum of {w_i c_i_x c_i_z} ",s)
        s = sum(w*cy*cy); call message(1,"TEST 3a: Sum of {w_i c_i_y c_i_y} ",s)
        s = sum(w*cy*cz); call message(1,"TEST 3a: Sum of {w_i c_i_y c_i_z} ",s)
        s = sum(w*cz*cz); call message(1,"TEST 3a: Sum of {w_i c_i_z c_i_z} ",s)

    end subroutine 
end module 
