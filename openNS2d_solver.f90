	subroutine Euler2d_solver(rho,u,v,p,xx,yy) !!
	
	include 'openNS3d.h'

	integer :: i,j,iter,ns

!	real(kind=OCFD_REAL_KIND) :: T

	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: u,v,p,rho
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: xx,yy
!	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) ::

    T=1.d0 
    ns=ceiling(T/dt)
	write(*,*) ns
	do iter=1,ns
		call RK4(rho,u,v,p)
		write(*,*) iter
		write(*,*) iter*dt

		if (iter==10000) then
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
			write(filename,"('Taylor',F5.1,'.plt')") iter*dt
			call write_data(rho,u,v,p,xx,yy,filename)
		endif

		if (mod(iter,25000)==0) then
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
			write(filename,"('Taylor',F5.1,'.plt')") iter*dt
			call write_data(rho,u,v,p,xx,yy,filename)
		endif
		
	end do




	end subroutine








!!!==========================================================================
!!!
!!!==========================================================================
subroutine computeR(rho,u,v,p,R_st,R_nd,R_rd,R_th)
	include 'openNS3d.h' 

!	real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:) :: Epos1,Epos2,Epos3,Epos4
!	real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:) :: Eneg1,Eneg2,Eneg3,Eneg4
!	real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:) :: Fpos1,Fpos2,Fpos3,Fpos4
!	real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:) :: Fneg1,Fneg2,Fneg3,Fneg4


!	allocate(Epos1(1-LAP:nx+LAP,1-LAP:ny+LAP),Epos2(1-LAP:nx+LAP,1-LAP:ny+LAP), &
!		Epos3(1-LAP:nx+LAP,1-LAP:ny+LAP),Epos4(1-LAP:nx+LAP,1-LAP:ny+LAP))
!	allocate(Eneg1(1-LAP:nx+LAP,1-LAP:ny+LAP),Eneg2(1-LAP:nx+LAP,1-LAP:ny+LAP), &
!		Eneg3(1-LAP:nx+LAP,1-LAP:ny+LAP),Eneg4(1-LAP:nx+LAP,1-LAP:ny+LAP))
!	allocate(Fpos1(1-LAP:nx+LAP,1-LAP:ny+LAP),Fpos2(1-LAP:nx+LAP,1-LAP:ny+LAP), &
!		Fpos3(1-LAP:nx+LAP,1-LAP:ny+LAP),Fpos4(1-LAP:nx+LAP,1-LAP:ny+LAP))
!	allocate(Fneg1(1-LAP:nx+LAP,1-LAP:ny+LAP),Fneg2(1-LAP:nx+LAP,1-LAP:ny+LAP), &
!		Fneg3(1-LAP:nx+LAP,1-LAP:ny+LAP),Fneg4(1-LAP:nx+LAP,1-LAP:ny+LAP))


	integer :: i,j

!	real(kind=OCFD_REAL_KIND) :: 

	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: u,v,p,rho
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: R_st,R_nd,R_rd,R_th
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Epos1,Epos2,Epos3,Epos4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Eneg1,Eneg2,Eneg3,Eneg4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Fpos1,Fpos2,Fpos3,Fpos4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Fneg1,Fneg2,Fneg3,Fneg4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: d2f_Epos1,d2f_Epos2,d2f_Epos3,d2f_Epos4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: d2f_Eneg1,d2f_Eneg2,d2f_Eneg3,d2f_Eneg4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: d2f_Fpos1,d2f_Fpos2,d2f_Fpos3,d2f_Fpos4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: d2f_Fneg1,d2f_Fneg2,d2f_Fneg3,d2f_Fneg4    
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: df_Epos1,df_Epos2,df_Epos3,df_Epos4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: df_Eneg1,df_Eneg2,df_Eneg3,df_Eneg4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: df_Fpos1,df_Fpos2,df_Fpos3,df_Fpos4
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: df_Fneg1,df_Fneg2,df_Fneg3,df_Fneg4 



	call VanLeerE(rho,u,v,p,Epos1,Epos2,Epos3,Epos4,& 
		Eneg1,Eneg2,Eneg3,Eneg4)
	call VanLeerF(rho,u,v,p,Fpos1,Fpos2,Fpos3,Fpos4,& 
		Fneg1,Fneg2,Fneg3,Fneg4)
	
	if (npx0 .ne. 1) then
		call check_x2d(Epos1)
		call check_x2d(Epos2)
		call check_x2d(Epos3)
		call check_x2d(Epos4)

		call check_x2d(Eneg1)
		call check_x2d(Eneg2)
		call check_x2d(Eneg3)
		call check_x2d(Eneg4)

		do j=1,ny
			call Du2Dx_PADE4(Epos1(:,j),d2f_Epos1(:,j),Iperiodic_X)
			call Du2Dx_PADE4(Eneg1(:,j),d2f_Eneg1(:,j),Iperiodic_X)
			call DuDx_UCC_UpWind(Epos1(:,j),d2f_Epos1(:,j),df_Epos1(:,j),Iperiodic_X)
			call DuDx_UCC_DownWind(Eneg1(:,j),d2f_Eneg1(:,j),df_Eneg1(:,j),Iperiodic_X)

			call Du2Dx_PADE4(Epos2(:,j),d2f_Epos2(:,j),Iperiodic_X)
			call Du2Dx_PADE4(Eneg2(:,j),d2f_Eneg2(:,j),Iperiodic_X)
			call DuDx_UCC_UpWind(Epos2(:,j),d2f_Epos2(:,j),df_Epos2(:,j),Iperiodic_X)
			call DuDx_UCC_DownWind(Eneg2(:,j),d2f_Eneg2(:,j),df_Eneg2(:,j),Iperiodic_X)

			call Du2Dx_PADE4(Epos3(:,j),d2f_Epos3(:,j),Iperiodic_X)
			call Du2Dx_PADE4(Eneg3(:,j),d2f_Eneg3(:,j),Iperiodic_X)
			call DuDx_UCC_UpWind(Epos3(:,j),d2f_Epos3(:,j),df_Epos3(:,j),Iperiodic_X)
			call DuDx_UCC_DownWind(Eneg3(:,j),d2f_Eneg3(:,j),df_Eneg3(:,j),Iperiodic_X)

			call Du2Dx_PADE4(Epos4(:,j),d2f_Epos4(:,j),Iperiodic_X)
			call Du2Dx_PADE4(Eneg4(:,j),d2f_Eneg4(:,j),Iperiodic_X)
			call DuDx_UCC_UpWind(Epos4(:,j),d2f_Epos4(:,j),df_Epos4(:,j),Iperiodic_X)
			call DuDx_UCC_DownWind(Eneg4(:,j),d2f_Eneg4(:,j),df_Eneg4(:,j),Iperiodic_X)
		end do

	elseif (npx0==1) then
		do j=1,ny
			call Du2Dx_PADE4_serial(Epos1(:,j),d2f_Epos1(:,j))
			call Du2Dx_PADE4_serial(Eneg1(:,j),d2f_Eneg1(:,j))
			call DuDx_UCC_UpWind_serial(Epos1(:,j),d2f_Epos1(:,j),df_Epos1(:,j))
			call DuDx_UCC_DownWind_serial(Eneg1(:,j),d2f_Eneg1(:,j),df_Eneg1(:,j))

			call Du2Dx_PADE4_serial(Epos2(:,j),d2f_Epos2(:,j))
			call Du2Dx_PADE4_serial(Eneg2(:,j),d2f_Eneg2(:,j))
			call DuDx_UCC_UpWind_serial(Epos2(:,j),d2f_Epos2(:,j),df_Epos2(:,j))
			call DuDx_UCC_DownWind_serial(Eneg2(:,j),d2f_Eneg2(:,j),df_Eneg2(:,j))

			call Du2Dx_PADE4_serial(Epos3(:,j),d2f_Epos3(:,j))
			call Du2Dx_PADE4_serial(Eneg3(:,j),d2f_Eneg3(:,j))
			call DuDx_UCC_UpWind_serial(Epos3(:,j),d2f_Epos3(:,j),df_Epos3(:,j))
			call DuDx_UCC_DownWind_serial(Eneg3(:,j),d2f_Eneg3(:,j),df_Eneg3(:,j))

			call Du2Dx_PADE4_serial(Epos4(:,j),d2f_Epos4(:,j))
			call Du2Dx_PADE4_serial(Eneg4(:,j),d2f_Eneg4(:,j))
			call DuDx_UCC_UpWind_serial(Epos4(:,j),d2f_Epos4(:,j),df_Epos4(:,j))
			call DuDx_UCC_DownWind_serial(Eneg4(:,j),d2f_Eneg4(:,j),df_Eneg4(:,j))
		end do
	endif

	if (npy0 .ne. 1) then
		call check_y2d(Fpos1)
		call check_y2d(Fpos2)
		call check_y2d(Fpos3)
		call check_y2d(Fpos4)

		call check_y2d(Fneg1)
		call check_y2d(Fneg2)
		call check_y2d(Fneg3)
		call check_y2d(Fneg4)

		do i=1,nx
			call Du2Dy_PADE4(Fpos1(i,:),d2f_Fpos1(i,:),Iperiodic_Y)
			call Du2Dy_PADE4(Fneg1(i,:),d2f_Fneg1(i,:),Iperiodic_Y)
			call DuDy_UCC_UpWind(Fpos1(i,:),d2f_Fpos1(i,:),df_Fpos1(i,:),Iperiodic_Y)
			call DuDy_UCC_DownWind(Fneg1(i,:),d2f_Fneg1(i,:),df_Fneg1(i,:),Iperiodic_Y)

			call Du2Dy_PADE4(Fpos2(i,:),d2f_Fpos2(i,:),Iperiodic_Y)
			call Du2Dy_PADE4(Fneg2(i,:),d2f_Fneg2(i,:),Iperiodic_Y)
			call DuDy_UCC_UpWind(Fpos2(i,:),d2f_Fpos2(i,:),df_Fpos2(i,:),Iperiodic_Y)
			call DuDy_UCC_DownWind(Fneg2(i,:),d2f_Fneg2(i,:),df_Fneg2(i,:),Iperiodic_Y)

			call Du2Dy_PADE4(Fpos3(i,:),d2f_Fpos3(i,:),Iperiodic_Y)
			call Du2Dy_PADE4(Fneg3(i,:),d2f_Fneg3(i,:),Iperiodic_Y)
			call DuDy_UCC_UpWind(Fpos3(i,:),d2f_Fpos3(i,:),df_Fpos3(i,:),Iperiodic_Y)
			call DuDy_UCC_DownWind(Fneg3(i,:),d2f_Fneg3(i,:),df_Fneg3(i,:),Iperiodic_Y)

			call Du2Dy_PADE4(Fpos4(i,:),d2f_Fpos4(i,:),Iperiodic_Y)
			call Du2Dy_PADE4(Fneg4(i,:),d2f_Fneg4(i,:),Iperiodic_Y)
			call DuDy_UCC_UpWind(Fpos4(i,:),d2f_Fpos4(i,:),df_Fpos4(i,:),Iperiodic_Y)
			call DuDy_UCC_DownWind(Fneg4(i,:),d2f_Fneg4(i,:),df_Fneg4(i,:),Iperiodic_Y)
		end do

	elseif (npy0==1) then
		do i=1,nx
			call Du2Dy_PADE4_serial(Fpos1(i,:),d2f_Fpos1(i,:))
			call Du2Dy_PADE4_serial(Fneg1(i,:),d2f_Fneg1(i,:))
			call DuDy_UCC_UpWind_serial(Fpos1(i,:),d2f_Fpos1(i,:),df_Fpos1(i,:))
			call DuDy_UCC_DownWind_serial(Fneg1(i,:),d2f_Fneg1(i,:),df_Fneg1(i,:))

			call Du2Dy_PADE4_serial(Fpos2(i,:),d2f_Fpos2(i,:))
			call Du2Dy_PADE4_serial(Fneg2(i,:),d2f_Fneg2(i,:))
			call DuDy_UCC_UpWind_serial(Fpos2(i,:),d2f_Fpos2(i,:),df_Fpos2(i,:))
			call DuDy_UCC_DownWind_serial(Fneg2(i,:),d2f_Fneg2(i,:),df_Fneg2(i,:))

			call Du2Dy_PADE4_serial(Fpos3(i,:),d2f_Fpos3(i,:))
			call Du2Dy_PADE4_serial(Fneg3(i,:),d2f_Fneg3(i,:))
			call DuDy_UCC_UpWind_serial(Fpos3(i,:),d2f_Fpos3(i,:),df_Fpos3(i,:))
			call DuDy_UCC_DownWind_serial(Fneg3(i,:),d2f_Fneg3(i,:),df_Fneg3(i,:))

			call Du2Dy_PADE4_serial(Fpos4(i,:),d2f_Fpos4(i,:))
			call Du2Dy_PADE4_serial(Fneg4(i,:),d2f_Fneg4(i,:))
			call DuDy_UCC_UpWind_serial(Fpos4(i,:),d2f_Fpos4(i,:),df_Fpos4(i,:))
			call DuDy_UCC_DownWind_serial(Fneg4(i,:),d2f_Fneg4(i,:),df_Fneg4(i,:))
		end do
	endif

	R_st=-(df_Epos1+df_Eneg1)-(df_Fpos1+df_Fneg1)
	R_nd=-(df_Epos2+df_Eneg2)-(df_Fpos2+df_Fneg2)
	R_rd=-(df_Epos3+df_Eneg3)-(df_Fpos3+df_Fneg3)
	R_th=-(df_Epos4+df_Eneg4)-(df_Fpos4+df_Fneg4)

	do j=1,ny
		do i=1,nx
			R_rd(i,j)=R_rd(i,j)+rho(i,j)
			R_th(i,j)=R_th(i,j)+rho(i,j)*v(i,j)
		end do
	end do
	

end subroutine




subroutine RK4(rho,u,v,p)
	include 'openNS3d.h'

	integer :: i,j

	real(kind=OCFD_REAL_KIND),parameter::gama=5.d0/3.d0

	real(kind=OCFD_REAL_KIND),parameter::M0=0.4d0,Mv=0.2d0,p0=1.d0/gama/M0**2.d0,rho0=1.d0,u0=1.d0

	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: u,v,p,rho

	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: R_st,R_nd,R_rd,R_th
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: R_st1,R_nd1,R_rd1,R_th1
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: R_st2,R_nd2,R_rd2,R_th2
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: R_st3,R_nd3,R_rd3,R_th3

	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Q_st,Q_nd,Q_rd,Q_th
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Q_st1,Q_nd1,Q_rd1,Q_th1
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Q_st2,Q_nd2,Q_rd2,Q_th2
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: Q_st3,Q_nd3,Q_rd3,Q_th3

!	gama=1.4d0

	call computeR(rho,u,v,p,R_st,R_nd,R_rd,R_th)


	do j=1,ny
		do i=1,nx
			Q_st(i,j)=rho(i,j)
			Q_nd(i,j)=rho(i,j)*u(i,j)
			Q_rd(i,j)=rho(i,j)*v(i,j)
			Q_th(i,j)=p(i,j)/(gama-1.d0)+rho(i,j)/2.d0*(u(i,j)**2.d0+v(i,j)**2.d0)
		end do
	end do

	Q_st1=Q_st+(dt/2.d0)*R_st
	Q_nd1=Q_nd+(dt/2.d0)*R_nd
	Q_rd1=Q_rd+(dt/2.d0)*R_rd
	Q_th1=Q_th+(dt/2.d0)*R_th

	do j=1,ny
		do i=1,nx
			rho(i,j)=Q_st1(i,j)
			u(i,j)=Q_nd1(i,j)/rho(i,j)
			v(i,j)=Q_rd1(i,j)/rho(i,j)
			p(i,j)=(gama-1.d0)*(Q_th1(i,j)-rho(i,j)/2.d0*(u(i,j)**2.d0+v(i,j)**2.d0))
		end do
	end do

!	if (npy==0) then
!		p(:,1)=1.d0
!		rho(:,1)=2.d0
!		u(:,1)=0.d0
!		v(:,1)=0.d0
!	endif

!	if (npy==npy0-1) then
!		p(:,ny)=2.5d0
!		rho(:,ny)=1.d0
!		u(:,ny)=0.d0
!		v(:,ny)=0.d0
!	endif


	call computeR(rho,u,v,p,R_st1,R_nd1,R_rd1,R_th1)


	Q_st2=Q_st+(dt/2.d0)*R_st1
	Q_nd2=Q_nd+(dt/2.d0)*R_nd1
	Q_rd2=Q_rd+(dt/2.d0)*R_rd1
	Q_th2=Q_th+(dt/2.d0)*R_th1

	do j=1,ny
		do i=1,nx
			rho(i,j)=Q_st2(i,j)
			u(i,j)=Q_nd2(i,j)/rho(i,j)
			v(i,j)=Q_rd2(i,j)/rho(i,j)
			p(i,j)=(gama-1.d0)*(Q_th2(i,j)-rho(i,j)/2.d0*(u(i,j)**2.d0+v(i,j)**2.d0))
		end do
	end do

!	if (npy==0) then
!		p(:,1)=1.d0
!		rho(:,1)=2.d0
!		u(:,1)=0.d0
!		v(:,1)=0.d0
!	endif

!	if (npy==npy0-1) then
!		p(:,ny)=2.5d0
!		rho(:,ny)=1.d0
!		u(:,ny)=0.d0
!		v(:,ny)=0.d0
!	endif

	call computeR(rho,u,v,p,R_st2,R_nd2,R_rd2,R_th2)

	Q_st3=Q_st+dt*R_st2
	Q_nd3=Q_nd+dt*R_nd2
	Q_rd3=Q_rd+dt*R_rd2
	Q_th3=Q_th+dt*R_th2
	
	do j=1,ny
		do i=1,nx
			rho(i,j)=Q_st3(i,j)
			u(i,j)=Q_nd3(i,j)/rho(i,j)
			v(i,j)=Q_rd3(i,j)/rho(i,j)
			p(i,j)=(gama-1.d0)*(Q_th3(i,j)-rho(i,j)/2.d0*(u(i,j)**2.d0+v(i,j)**2.d0))
		end do
	end do
	
!	if (npy==0) then
!		p(:,1)=1.d0
!		rho(:,1)=2.d0
!		u(:,1)=0.d0
!		v(:,1)=0.d0
!	endif

!	if (npy==npy0-1) then
!		p(:,ny)=2.5d0
!		rho(:,ny)=1.d0
!		u(:,ny)=0.d0
!		v(:,ny)=0.d0
!	endif


	call computeR(rho,u,v,p,R_st3,R_nd3,R_rd3,R_th3)

	Q_st=Q_st+(dt/6.d0)*(R_st+2.d0*R_st1+2.d0*R_st2+R_st3)
	Q_nd=Q_nd+(dt/6.d0)*(R_nd+2.d0*R_nd1+2.d0*R_nd2+R_nd3)
	Q_rd=Q_rd+(dt/6.d0)*(R_rd+2.d0*R_rd1+2.d0*R_rd2+R_rd3)
	Q_th=Q_th+(dt/6.d0)*(R_th+2.d0*R_th1+2.d0*R_th2+R_th3)
	
	do j=1,ny
		do i=1,nx
			rho(i,j)=Q_st(i,j)
			u(i,j)=Q_nd(i,j)/rho(i,j)
			v(i,j)=Q_rd(i,j)/rho(i,j)
			p(i,j)=(gama-1.d0)*(Q_th(i,j)-rho(i,j)/2.d0*(u(i,j)**2.d0+v(i,j)**2.d0))
		end do
	end do

!	if (npy==0) then
!		p(:,1)=1.d0
!		rho(:,1)=2.d0
!		u(:,1)=0.d0
!		v(:,1)=0.d0
!	endif
!
!	if (npy==npy0-1) then
!		p(:,ny)=2.5d0
!		rho(:,ny)=1.d0
!		u(:,ny)=0.d0
!		v(:,ny)=0.d0
!	endif



end subroutine RK4