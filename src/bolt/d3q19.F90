module d3q19mod
    use kind_parameters, only: rkind, clen 
    use exits, only: GracefulExit, message, check_exit
    use decomp_2d
    use constants, only: one, zero, two 
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

    type :: d3q19
        private
        
        type(decomp_info) :: gp
        integer :: nx, ny, nz

        character(len=clen) :: inputdir, outputdir

        real(rkind), dimension(:,:,:,:), allocatable :: f, feq

        real(rkind) :: delta_t, delta_x, delta_u
        real(rkind) :: Re

        real(rkind), dimension(:,:,:), allocatable :: rho, ux, uy, uz, nu, tau, Fx, Fy, Fz

        real(rkind), dimension(3,3,nvels) :: QTensor
        real(rkind), dimension(:,:,:), allocatable :: PiTensor

        contains
            procedure :: init
            procedure :: destroy

            procedure :: collide 
            procedure :: stream

            procedure :: compute_macroscopic 

    end type


    logical, parameter :: RestrictSingleCore = .true. 
contains

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


    subroutine compute_feq(this)
        class(d3q19), intent(inout) :: this 
        
        integer :: i, j, k, idx 
        real(rkind) :: first, second


        do idx = 1,nvels
            do k = 1,this%gp%xsz(3)
                do j = 1,this%gp%xsz(2)
                    !$omp simd
                    do i = 1,this%gp%xsz(1)
                        first = cx(idx)*this%ux(i,j,k) + cy(idx)*this%uy(i,j,k) + cz(idx)*this%uz(i,j,k)

                        second =     this%ux(i,j,k)*this%ux(i,j,k)*this%Qtensor(1,1,idx) &
                             & +     this%uy(i,j,k)*this%uy(i,j,k)*this%Qtensor(2,2,idx) &
                             & +     this%uz(i,j,k)*this%uz(i,j,k)*this%Qtensor(3,3,idx) &
                             & + two*this%ux(i,j,k)*this%uy(i,j,k)*this%Qtensor(1,2,idx) &
                             & + two*this%ux(i,j,k)*this%uz(i,j,k)*this%Qtensor(1,3,idx) &
                             & + two*this%uy(i,j,k)*this%uz(i,j,k)*this%Qtensor(2,3,idx) 
                    
                        this%feq(i,j,k,idx) = w(idx)*this%rho(i,j,k)*(one + onebycsq*first + oneby2c4*second)

                    end do 
                end do 
            end do 
        end do

    end subroutine 


    subroutine compute_macroscopic(this)
        class(d3q19), intent(inout) :: this
        integer :: i, j, k, onebyrho

        do k = 1,this%gp%xsz(3)
            do j = 1,this%gp%xsz(2)
                !$omp simd
                do i = 1,this%gp%xsz(1)
                    this%rho(i,j,k) = this%f(i,j,k,1)  + this%f(i,j,k,2)  + this%f(i,j,k,3)  + this%f(i,j,k,4)  &
                                 &  + this%f(i,j,k,5)  + this%f(i,j,k,6)  + this%f(i,j,k,7)  + this%f(i,j,k,8)  &
                                 &  + this%f(i,j,k,9)  + this%f(i,j,k,10) + this%f(i,j,k,11) + this%f(i,j,k,12) &
                                 &  + this%f(i,j,k,13) + this%f(i,j,k,14) + this%f(i,j,k,15) + this%f(i,j,k,16) &
                                 &  + this%f(i,j,k,17) + this%f(i,j,k,18) + this%f(i,j,k,19)

                    onebyrho = one/this%rho(i,j,k)  

                    this%ux(i,j,k)  = (this%f(i,j,k,1)  - this%f(i,j,k,2)  + this%f(i,j,k,7)  - this%f(i,j,k,8)  &
                                 &  +  this%f(i,j,k,9)  - this%f(i,j,k,10) + this%f(i,j,k,13) - this%f(i,j,k,14) &
                                 &  +  this%f(i,j,k,15) - this%f(i,j,k,16))*onebyrho  

                    this%uy(i,j,k)  = (this%f(i,j,k,3)  - this%f(i,j,k,4)  + this%f(i,j,k,7)  - this%f(i,j,k,8)  &
                                 &  +  this%f(i,j,k,11) - this%f(i,j,k,12) - this%f(i,j,k,13) + this%f(i,j,k,14) &
                                 &  +  this%f(i,j,k,17) - this%f(i,j,k,18))*onebyrho  
                    
                    this%uz(i,j,k)  = (this%f(i,j,k,5)  - this%f(i,j,k,6)  + this%f(i,j,k,9)  - this%f(i,j,k,10) &
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
