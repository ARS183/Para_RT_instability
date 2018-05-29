
program shielded_vortex

	implicit none
	include 'openNS3d.h'	

	real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:) :: rho,u,v,p,xx,yy
	real(kind=OCFD_REAL_KIND) :: time,time_start,time_end

	real(kind=OCFD_REAL_KIND),parameter::xo=4.d0,yo=0.d0
	real(kind=OCFD_REAL_KIND),parameter::M0=0.4d0,Mv=0.2d0,gama=1.4d0,p0=1.d0/gama/M0**2.d0,rho0=1.d0,u0=1.d0
	real(kind=OCFD_REAL_KIND)::Ms=M0,ps=p0

	call MPI_INIT(ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,np_size,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)

	

	call read_parameter()
	
	call part2d()	
	allocate(u(1-LAP:nx+LAP,1-LAP:ny+LAP),v(1-LAP:nx+LAP,1-LAP:ny+LAP),p(1-LAP:nx+LAP,1-LAP:ny+LAP), &
		xx(1-LAP:nx+LAP,1-LAP:ny+LAP),yy(1-LAP:nx+LAP,1-LAP:ny+LAP),rho(1-LAP:nx+LAP,1-LAP:ny+LAP))

	time=0.d0
	call init_RT(rho,u,v,p,xx,yy) !!
	!filename='field_init.plt'
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	
!	time_start=MPI_WTIME()
	
	write(filename,"('Taylor',F4.1,'.plt')") time
	call write_data(rho,u,v,p,xx,yy,filename)

	
	time_start=MPI_WTIME()
	call Euler2d_solver(rho,u,v,p,xx,yy)
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	time_end=MPI_WTIME()

	if(my_id.eq.0) then
		open(10,file='walltime.dat',status='unknown')	
 		write(10,"('Time=',E16.6)")time_end-time_start
	endif

	deallocate(u,v,p,xx,yy,rho)  !!

	call MPI_FINALIZE(ierr)

end program shielded_vortex


	subroutine read_parameter()
	include 'openNS3d.h'

	Nparameters=0
	Rparameters=0.d0	

	if (my_id .eq. 0) then
	open(100,file='openNS2d.in')
	
	read(100,*)
	read(100,*)
	read(100,*)nx_global,ny_global
	read(100,*)
	read(100,*)npx0,npy0
	read(100,*)
	read(100,*)slx,sly
	read(100,*)
	read(100,*)Iperiodic_X,Iperiodic_Y
	read(100,*)
	read(100,*)Re,dt,rf,eps
	read(100,*)
	read(100,*)Istep_show,Istep_save
	read(100,*)
	read(100,*)NUM_METHOD_CONV,NUM_METHOD_OTH
	

	close(100)
	
	endif

	Nparameters(1)=nx_global
	Nparameters(2)=ny_global
	Nparameters(3)=npx0
	Nparameters(4)=npy0
	Nparameters(5)=Istep_show
	Nparameters(6)=Istep_save
	Nparameters(7)=Iperiodic_X
	Nparameters(8)=Iperiodic_Y
	Nparameters(9)=NUM_METHOD_CONV
	Nparameters(10)=NUM_METHOD_OTH

	Rparameters(1)=slx
	Rparameters(2)=sly
	Rparameters(3)=Re
	Rparameters(4)=dt
	Rparameters(5)=rf
	Rparameters(6)=eps

	call MPI_BCAST(Nparameters,NparaMax,MPI_integer,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(Rparameters,RparaMax,OCFD_DATA_TYPE,0,MPI_COMM_WORLD,ierr)
	
	nx_global=Nparameters(1)
	ny_global=Nparameters(2)
	npx0=Nparameters(3)
	npy0=Nparameters(4)
	Istep_show=Nparameters(5)
	Istep_save=Nparameters(6)
	Iperiodic_X=Nparameters(7)
	Iperiodic_Y=Nparameters(8)
	NUM_METHOD_CONV=Nparameters(9)
	NUM_METHOD_OTH=Nparameters(10)


	
	slx=Rparameters(1)
	sly=Rparameters(2)
	Re=Rparameters(3)
	dt=Rparameters(4)
	rf=Rparameters(5)
	eps=Rparameters(6)

	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	

	end subroutine read_parameter
