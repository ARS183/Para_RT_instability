!	implicit none
	include 'mpif.h'

!===double precision computation or single precision computation
    integer,parameter::OCFD_REAL_KIND=8, OCFD_DATA_TYPE=MPI_DOUBLE_PRECISION   !double precision computation 
!   integer,parameter::OCFD_REAL_KIND=4, OCFD_DATA_TYPE=MPI_REAL               !single presision  computation

!===program parameter
	integer,parameter :: npMax=1024  ! the maximun number of computation nodes permitted
	integer,parameter :: NparaMax=100          ! the maximun number of integer input parameter
	integer,parameter :: RparaMax=100          ! the maximun number of real input parameter
	real(kind=OCFD_REAL_KIND) :: eps,eps_stable
	character*100 filename

	common /progpara/eps,eps_stable

!===partition===
	integer :: nx_global,ny_global,nz_global,npx0,npy0,npz0,npx,npy,npz,nx,ny,nz,nx1,nx2,ny1,ny2,nz1,nz2
	integer :: ID_XM1,ID_XP1,ID_YM1,ID_YP1,ID_ZM1,ID_ZP1,MPI_COMM_X,MPI_COMM_Y,MPI_COMM_Z,Iperiodic_X,Iperiodic_Y,Iperiodic_Z, &
	    MPI_COMM_XY,MPI_COMM_XZ,MPI_COMM_YZ
	integer :: TYPE_LAPX1,TYPE_LAPY1,TYPE_LAPZ1,TYPE_LAPX2,TYPE_LAPY2,TYPE_LAPZ2
	real(kind=OCFD_REAL_KIND) :: slx,sly,slz,hx,hy,hz
	integer,parameter :: OCFD_Barrier_level=1
	integer,parameter :: MSG_BLOCK_SIZE=0
	! 以上两个参数意义见mpi，其他模块可能暂不可用
	integer,dimension(0:npMax-1) :: i_offset,j_offset,k_offset,i_nn,j_nn,k_nn
	
	common /partition/i_offset,j_offset,k_offset,i_nn,j_nn,k_nn,npx,npy,npz,nx,ny,nz,nx1,nx2,ny1,ny2,nz1,nz2,  &
			ID_XM1,ID_XP1,ID_YM1,ID_YP1,ID_ZM1,ID_ZP1,MPI_COMM_X,MPI_COMM_Y,MPI_COMM_Z,MPI_COMM_XY,MPI_COMM_XZ,MPI_COMM_YZ,  &
			nx_global,ny_global,nz_global,npx0,npy0,npz0,Iperiodic_X,Iperiodic_Y,Iperiodic_Z, &
			TYPE_LAPX1,TYPE_LAPY1,TYPE_LAPZ1,TYPE_LAPX2,TYPE_LAPY2,TYPE_LAPZ2
	common /partition1/slx,sly,slz,hx,hy,hz

!===problem variables

!===problem parameter===

	real(kind=OCFD_REAL_KIND) :: Re,dt,rf,T,Rem,Ha
	integer :: Istep_show,Istep_save
	integer,dimension(1:NparaMax) :: Nparameters
	real(kind=OCFD_REAL_KIND),dimension(1:RparaMax) :: Rparameters
	integer,parameter :: LAP=4

	common /propara/Re,dt,rf,T,Rem,Ha,Istep_show,Istep_save,Nparameters,Rparameters

!===numerical method===
	integer :: NUM_METHOD_CONV,NUM_METHOD_OTH,NUM_TIME,i_case,i_grid
	common /numerical/ NUM_METHOD_CONV,NUM_METHOD_OTH,NUM_TIME,i_case,i_grid
	integer,parameter:: OCFD_NUMERICAL_AECDS_2nd=1,OCFD_NUMERICAL_AECDS_4th=2,	&
	OCFD_NUMERICAL_AECDS_6th=3,OCFD_NUMERICAL_AECDS_o2nd=4,OCFD_NUMERICAL_AECDS_o4th=5,	&
	OCFD_NUMERICAL_AECDS_o6th=6,OCFD_NUMERICAL_WENO_5th=7,OCFD_NUMERICAL_WENO_7th=8,	&
	OCFD_POSSION_2nd=1,OCFD_POSSION_4th=2,OCFD_POSSION_6th=3,DNS=1,LES=2


!===MPI variables===
	integer :: ierr,np_size,my_id
	integer :: STATUS(MPI_STATUS_SIZE),request1,request2

	common /mpivar/np_size,my_id
