!===BCs for different problems
	subroutine OCFD_DFX_BOUND_CHECK_2d(u,f,NUM_METHOD)
	include 'openNS3d.h'
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: u,f
	integer i,j,NUM_METHOD

    if (Iperiodic_X .ne. 1) then !!!!

	if (npx .eq. 0) then
	do j=1,ny
	call OCFD_DF_BOUND(u(:,j),f(:,j),nx,hx,NUM_METHOD,1)
	enddo
	endif

	if (npx .eq. npx0-1) then
	do j=1,ny
	call OCFD_DF_BOUND(u(:,j),f(:,j),nx,hx,NUM_METHOD,2)
	enddo
	endif

	endif


	end subroutine

!============================
	subroutine OCFD_DFY_BOUND_CHECK_2d(u,f,NUM_METHOD)
	include 'openNS3d.h'
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: u,f
	integer i,j,NUM_METHOD

        
	if (Iperiodic_Y .ne. 1) then !!!
	if (npy .eq. 0) then
	do i=1,nx
	call OCFD_DF_BOUND(u(i,:),f(i,:),ny,hy,NUM_METHOD,1)
	enddo
	endif

	if (npy .eq. npy0-1) then
	do i=1,nx
	call OCFD_DF_BOUND(u(i,:),f(i,:),ny,hy,NUM_METHOD,2)
	enddo
	endif
	
	endif

	end subroutine

	subroutine OCFD_DFX_BOUND_CHECK(u,f,NUM_METHOD)
	include 'openNS3d.h'
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP) :: u,f
	integer i,j,k,NUM_METHOD

    if (Iperiodic_X .ne. 1) then !!!!

	if (npx .eq. 0) then
	do k=1,nz
	do j=1,ny
	call OCFD_DF_BOUND(u(:,j,k),f(:,j,k),nx,hx,NUM_METHOD,1)
	enddo
	enddo
	endif

	if (npx .eq. npx0-1) then
	do k=1,nz
	do j=1,ny
	call OCFD_DF_BOUND(u(:,j,k),f(:,j,k),nx,hx,NUM_METHOD,2)
	enddo
	enddo
	endif

	endif


	end subroutine

!============================
	subroutine OCFD_DFY_BOUND_CHECK(u,f,NUM_METHOD)
	include 'openNS3d.h'
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP) :: u,f
	integer i,j,k,NUM_METHOD

        
	if (Iperiodic_Y .ne. 1) then !!!
	if (npy .eq. 0) then
	do i=1,nx
	do k=1,nz
	call OCFD_DF_BOUND(u(i,:,k),f(i,:,k),ny,hy,NUM_METHOD,1)
	enddo
	enddo
	endif

	if (npy .eq. npy0-1) then
	do i=1,nx
	do k=1,nz
	call OCFD_DF_BOUND(u(i,:,k),f(i,:,k),ny,hy,NUM_METHOD,2)
	enddo
	enddo
	endif
	
	endif

	end subroutine

!============================
	subroutine OCFD_DFZ_BOUND_CHECK(u,f,NUM_METHOD)
	include 'openNS3d.h'
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP) :: u,f
	integer i,j,k,NUM_METHOD

        
	if (Iperiodic_Z .ne. 1) then !!!
	if (npz .eq. 0) then
	do i=1,nx
	do j=1,ny
	call OCFD_DF_BOUND(u(i,j,:),f(i,j,:),nz,hz,NUM_METHOD,1)
	enddo
	enddo
	endif

	if (npz .eq. npz0-1) then
	do i=1,nx
	do j=1,ny
	call OCFD_DF_BOUND(u(i,j,:),f(i,j,:),nz,hz,NUM_METHOD,2)
	enddo
	enddo
	endif
	
	endif

	end subroutine

!============================

	subroutine OCFD_DF_BOUND(u,f,n,h,NUM_METHOD,flag)
	include 'openNS3d.h'
	real(kind=OCFD_REAL_KIND)::h
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)
	integer n,NUM_METHOD,flag

	if (NUM_METHOD .eq. OCFD_NUMERICAL_AECDS_2nd) then
	call OCFD_DF_BOUND_2nd(u,f,n,h,flag)
	goto 100
	
	else if(NUM_METHOD .eq. OCFD_NUMERICAL_AECDS_4th .or. NUM_METHOD .eq. OCFD_NUMERICAL_AECDS_o2nd) then
	call OCFD_DF_BOUND_4th(u,f,n,h,flag)
	goto 100

	else if(NUM_METHOD .eq. OCFD_NUMERICAL_AECDS_6th .or. NUM_METHOD .eq. OCFD_NUMERICAL_AECDS_o4th) then
	call OCFD_DF_BOUND_6th(u,f,n,h,flag)
	goto 100


	else
	      print*, 'This Numerical Method is not supported'
		  stop     
	endif


100	continue


	end subroutine

!=================================


	subroutine OCFD_DF2X_BOUND_CHECK(u,f,NUM_METHOD)
	include 'openNS3d.h'
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP) :: u,f
	integer i,j,k,NUM_METHOD

        if (Iperiodic_X .ne. 1) then !!!!

	if (npx .eq. 0) then
	do j=1,ny
	do k=1,nz
	call OCFD_DF2_BOUND(u(:,j,k),f(:,j,k),nx,hx,NUM_METHOD,1)
	enddo
	enddo
	endif

	if (npx .eq. npx0-1) then
	do j=1,ny
	do k=1,nz
	call OCFD_DF2_BOUND(u(:,j,k),f(:,j,k),nx,hx,NUM_METHOD,2)
	enddo
	enddo
	endif

	endif


	end subroutine

!============================
	subroutine OCFD_DF2Y_BOUND_CHECK_2d(u,f,NUM_METHOD)
	include 'openNS3d.h'
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: u,f
	integer i,j,k,NUM_METHOD

        
	if (Iperiodic_Y .ne. 1) then !!!
	if (npy .eq. 0) then
	do i=1,nx
	call OCFD_DF2_BOUND(u(i,:),f(i,:),ny,hy,NUM_METHOD,1)
	enddo
	endif

	if (npy .eq. npy0-1) then
	do i=1,nx
	call OCFD_DF2_BOUND(u(i,:),f(i,:),ny,hy,NUM_METHOD,2)
	enddo
	endif
	
	endif

	end subroutine


!=================================


	subroutine OCFD_DF2X_BOUND_CHECK_2d(u,f,NUM_METHOD)
	include 'openNS3d.h'
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: u,f
	integer i,j,NUM_METHOD

        if (Iperiodic_X .ne. 1) then !!!!

	if (npx .eq. 0) then
	do j=1,ny
	call OCFD_DF2_BOUND(u(:,j),f(:,j),nx,hx,NUM_METHOD,1)
	enddo
	endif

	if (npx .eq. npx0-1) then
	do j=1,ny
	call OCFD_DF2_BOUND(u(:,j),f(:,j),nx,hx,NUM_METHOD,2)
	enddo
	endif

	endif


	end subroutine

!============================
	subroutine OCFD_DF2Y_BOUND_CHECK(u,f,NUM_METHOD)
	include 'openNS3d.h'
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP) :: u,f
	integer i,j,k,NUM_METHOD

        
	if (Iperiodic_Y .ne. 1) then !!!
	if (npy .eq. 0) then
	do i=1,nx
	do k=1,nz
	call OCFD_DF2_BOUND(u(i,:,k),f(i,:,k),ny,hy,NUM_METHOD,1)
	enddo
	enddo
	endif

	if (npy .eq. npy0-1) then
	do i=1,nx
	do k=1,nz
	call OCFD_DF2_BOUND(u(i,:,k),f(i,:,k),ny,hy,NUM_METHOD,2)
	enddo
	enddo
	endif
	
	endif

	end subroutine

!============================
	subroutine OCFD_DF2Z_BOUND_CHECK(u,f,NUM_METHOD)
	include 'openNS3d.h'
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP,1-LAP:nz+LAP) :: u,f
	integer i,j,k,NUM_METHOD

        
	if (Iperiodic_Z .ne. 1) then !!!
	if (npz .eq. 0) then
	do i=1,nx
	do j=1,ny
	call OCFD_DF2_BOUND(u(i,j,:),f(i,j,:),nz,hz,NUM_METHOD,1)
	enddo
	enddo
	endif

	if (npz .eq. npz0-1) then
	do i=1,nx
	do j=1,ny
	call OCFD_DF2_BOUND(u(i,j,:),f(i,j,:),nz,hz,NUM_METHOD,2)
	enddo
	enddo
	endif
	
	endif

	end subroutine

!============================

	subroutine OCFD_DF2_BOUND(u,f,n,h,NUM_METHOD,flag)
	include 'openNS3d.h'
	real(kind=OCFD_REAL_KIND)::h
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)
	integer n,NUM_METHOD,flag

	if (NUM_METHOD .eq. OCFD_NUMERICAL_AECDS_2nd) then
	call OCFD_DF2_BOUND_2nd(u,f,n,h,flag)
	goto 100
	
	else if(NUM_METHOD .eq. OCFD_NUMERICAL_AECDS_4th) then
	call OCFD_DF2_BOUND_4th(u,f,n,h,flag)
	goto 100

	else if(NUM_METHOD .eq. OCFD_NUMERICAL_AECDS_6th) then
	call OCFD_DF2_BOUND_6th(u,f,n,h,flag)
	goto 100


	else
	      print*, 'This Numerical Method is not supported'
		  stop     
	endif


100	continue


	end subroutine

!=================================
	subroutine OCFD_DF_BOUND_2nd(u,f,n,h,flag)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)
	integer n,flag

	if (flag.eq.1) then
	f(1)=(-u(3)+4.d0*u(2)-3.d0*u(1))/2.d0/h
	elseif (flag.eq.2) then
	f(n)=-(-u(n-2)+4.d0*u(n-1)-3.d0*u(n))/2.d0/h
	endif


	end subroutine

!=================================
	subroutine OCFD_DF_BOUND_4th(u,f,n,h,flag)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)
	integer n,flag

	if (flag.eq.1) then
	f(1)=(-3.d0*u(5)+16.d0*u(4)-36.d0*u(3)+48.d0*u(2)-25.d0*u(1))/12.d0/h
	f(2)=(1.d0*u(5)-6.d0*u(4)+18.d0*u(3)-10.d0*u(2)-3.d0*u(1))/12.d0/h
	elseif (flag.eq.2) then
	f(n)=-(-3.d0*u(n-4)+16.d0*u(n-3)-36.d0*u(n-2)+48.d0*u(n-1)-25.d0*u(n))/12.d0/h
	f(n-1)=-(1.d0*u(n-4)-6.d0*u(n-3)+18.d0*u(n-2)-10.d0*u(n-1)-3.d0*u(n))/12.d0/h
	endif


	end subroutine

!=================================
	subroutine OCFD_DF_BOUND_6th_1(u,f,n,h,flag)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)
	integer n,flag

	if (flag.eq.1) then
	f(1)=(-3.d0*u(5)+16.d0*u(4)-36.d0*u(3)+48.d0*u(2)-25.d0*u(1))/12.d0/h
	f(2)=(1.d0*u(5)-6.d0*u(4)+18.d0*u(3)-10.d0*u(2)-3.d0*u(1))/12.d0/h
	f(3)=(-1.d0*u(5)+8.d0*u(4)-8.d0*u(2)+1.d0*u(1))/12.d0/h
	elseif (flag.eq.2) then
	f(n)=-(-3.d0*u(n-4)+16.d0*u(n-3)-36.d0*u(n-2)+48.d0*u(n-1)-25.d0*u(n))/12.d0/h
	f(n-1)=-(1.d0*u(n-4)-6.d0*u(n-3)+18.d0*u(n-2)-10.d0*u(n-1)-3.d0*u(n))/12.d0/h
	f(n-2)=-(-1.d0*u(n-4)+8.d0*u(n-3)-8.d0*u(n-1)+1.d0*u(n))/12.d0/h
	endif


	end subroutine

!==============================
	subroutine OCFD_DF_BOUND_6th(u,f,n,h,flag)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)
	integer n,flag

	if (flag.eq.1) then
	f(1)=(-147.d0*u(1)+360.d0*u(2)-450.d0*u(3)+400.d0*u(4)-225.d0*u(5)+72.d0*u(6)-10.d0*u(7))/60.d0/h
	f(2)=(-10.d0*u(1)-77.d0*u(2)+150.d0*u(3)-100.d0*u(4)+50.d0*u(5)-15.d0*u(6)+2.d0*u(7))/60.d0/h
	f(3)=(2.d0*u(1)-24.d0*u(2)-35.d0*u(3)+80.d0*u(4)-30.d0*u(5)+8.d0*u(6)-1.d0*u(7))/60.d0/h
	elseif (flag.eq.2) then
	f(n)=-(-147.d0*u(n)+360.d0*u(n-1)-450.d0*u(n-2)+400.d0*u(n-3)-225.d0*u(n-4)+72.d0*u(n-5)-10.d0*u(n-6))/60.d0/h
	f(n-1)=-(-10.d0*u(n)-77.d0*u(n-1)+150.d0*u(n-2)-100.d0*u(n-3)+50.d0*u(n-4)-15.d0*u(n-5)+2.d0*u(n-6))/60.d0/h
	f(n-2)=-(2.d0*u(n)-24.d0*u(n-1)-35.d0*u(n-2)+80.d0*u(n-3)-30.d0*u(n-4)+8.d0*u(n-5)-1.d0*u(n-6))/60.d0/h
	endif


	end subroutine

!=================================
	subroutine OCFD_DF2_BOUND_2nd(u,f,n,h,flag)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)
	integer n,flag

	if (flag.eq.1) then
	f(1)=(-u(4)+4.d0*u(3)-5.d0*u(2)+2.d0*u(1))/h/h
	elseif (flag.eq.2) then
	f(n)=(-u(n-3)+4.d0*u(n-2)-5.d0*u(n-1)+2.d0*u(n))/h/h
	endif


	end subroutine

!=================================
	subroutine OCFD_DF2_BOUND_4th(u,f,n,h,flag)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)
	integer n,flag

!	if (flag.eq.1) then
!	f(1)=(-10.d0*u(6)+61.d0*u(5)-156.d0*u(4)+214.d0*u(3)-154.d0*u(2)+45.d0*u(1))/12.d0/h/h
!	f(2)=(u(6)-6.d0*u(5)+14.d0*u(4)-4.d0*u(3)-15.d0*u(2)+10.d0*u(1))/12.d0/h/h
!	elseif (flag.eq.2) then
!	f(n)=(-10.d0*u(n-5)+61.d0*u(n-4)-156.d0*u(n-3)+214.d0*u(n-2) &
!		-154.d0*u(n-1)+45.d0*u(n))/12.d0/h/h
!	f(n-1)=(u(n-5)-6.d0*u(n-4)+14.d0*u(n-3)-4.d0*u(n-2)-15.d0*u(n-1)+10.d0*u(n))/12.d0/h/h
!	endif


	if (flag.eq.1) then
	f(1)=(-10.d0*u(6)+61.d0*u(5)-156.d0*u(4)+214.d0*u(3)-154.d0*u(2)+45.d0*u(1))/12.d0/h/h
	f(2)=(u(6)-6.d0*u(5)+14.d0*u(4)-4.d0*u(3)-15.d0*u(2)+10.d0*u(1))/12.d0/h/h
	f(3)=(-u(1)+16.d0*u(2)-30.d0*u(3)+16.d0*u(4)-u(5))/12.d0/h/h
	elseif (flag.eq.2) then
	f(n)=(-10.d0*u(n-5)+61.d0*u(n-4)-156.d0*u(n-3)+214.d0*u(n-2) &
		-154.d0*u(n-1)+45.d0*u(n))/12.d0/h/h
	f(n-1)=(u(n-5)-6.d0*u(n-4)+14.d0*u(n-3)-4.d0*u(n-2)-15.d0*u(n-1)+10.d0*u(n))/12.d0/h/h
    f(n-2)=(-u(n)+16.d0*u(n-1)-30.d0*u(n-2)+16.d0*u(n-3)-u(n-4))/12.d0/h/h
	endif


	end subroutine

!==============================
	subroutine OCFD_DF2_BOUND_6th(u,f,n,h,flag)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)
	integer n,flag

	if (flag.eq.1) then
	f(1)=(-3.d0*u(5)+16.d0*u(4)-36.d0*u(3)+48.d0*u(2)-25.d0*u(1))/12.d0/h
	f(2)=(1.d0*u(5)-6.d0*u(4)+18.d0*u(3)-10.d0*u(2)-3.d0*u(1))/12.d0/h
	elseif (flag.eq.2) then
	f(n)=-(-3.d0*u(n-4)+16.d0*u(n-3)-36.d0*u(n-2)+48.d0*u(n-1)-25.d0*u(n))/12.d0/h
	f(n-1)=-(1.d0*u(n-4)-6.d0*u(n-3)+18.d0*u(n-2)-10.d0*u(n-1)-3.d0*u(n))/12.d0/h
	endif


	end subroutine

	

!--by³Â½ðÇ¿ 2018´º
!!!==========================DF_BOUND_UCC45==========================!!!
!! Boundary
subroutine OCFD_DF_BOUND_UCC45(u,f,n,h,flag)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h,pi
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)
	integer n,flag

	pi=4.d0*datan(1.d0)

	if (flag.eq.1) then
!		f(1)=-2.d0*pi*(dsin(2.d0*pi*0.d0)-dcos(2.d0*pi*0.d0))
!		f(2)=-2.d0*pi*(dsin(2.d0*pi*h)-dcos(2.d0*pi*h))
!		f(1)=-2.d0*(0.d0-5.d0)*exp(-(0.d0-5.d0)**2)
!		f(2)=-2.d0*(h-5.d0)*exp(-(h-5.d0)**2)
		f(1)=0.d0
		f(2)=0.d0
	elseif (flag.eq.2) then
!		f(n)=-2.d0*pi*(dsin(2.d0*pi*1.d0)-dcos(2.d0*pi*1.d0))
!		f(n-1)=-2.d0*pi*(dsin(2.d0*pi*(1.d0-h))-dcos(2.d0*pi*(1.d0-h)))
!		f(n)=-2.d0*(10.d0-5.d0)*exp(-(10.d0-5.d0)**2)
!		f(n-1)=-2.d0*((10.d0-h)-5.d0)*exp(-((10.d0-h)-5.d0)**2)
		f(n)=0.d0
		f(n-1)=0.d0
	endif
end subroutine



!! Subboundary
subroutine OCFD_DF_SB_UCC45_P(u,f,n,h)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h
    integer n
	real(kind=OCFD_REAL_KIND)::a(-3:2)
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)

	
	a(-3)=-0.0314739841972651d0
	a(-2)=0.2407032543196588d0
	a(-1)=-0.9814065086393177d0
	a(0)=0.3147398419726509d0
	a(1)=0.5092967456803412d0
	a(2)=-0.051859349136068184d0

	f(1)=(a(-3)*u(-2)+a(-2)*u(-1)+a(-1)*u(0) &
	+a(0)*u(1)+a(1)*u(2)+a(2)*u(3))/h
	
end subroutine


subroutine OCFD_DF_SB_UCC45_M(u,f,n,h)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h
    integer n
	real(kind=OCFD_REAL_KIND)::a(-2:3)
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)

	
	a(3)=0.0314739841972651d0
	a(2)=-0.2407032543196588d0
	a(1)=0.9814065086393177d0
	a(0)=-0.3147398419726509d0
	a(-1)=-0.5092967456803412d0
	a(-2)=0.051859349136068184d0

	f(n)=(a(-2)*u(n-2)+a(-1)*u(n-1) &
	+a(0)*u(n)+a(1)*u(n+1)+a(2)*u(n+2)+a(3)*u(n+3))/h
	
end subroutine

!!!===============================================================!!!


!!!================================================================!!!
!DF_SB2_UCC45,subboundary,cal by the opposite scheme                 !
!!!================================================================!!!
subroutine OCFD_DF_SB2_UCC45_P(u,f,s,n,h)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h,sig,esi,h2,h4,h5,h6
    integer n
	real(kind=OCFD_REAL_KIND)::r(1:2),d6f(1:2)
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)
	real(kind=OCFD_REAL_KIND)::s(1-LAP:n+LAP)

	sig=0.54395d0
!	esi=-0.01385977393639451d0
	esi=-0.01820010801241671d0
	h2=h*h
	h4=h2*h2
	h5=h4*h
	h6=h4*h2

	d6f(1)=(u(-2)-6.d0*u(-1)+15.d0*u(0)-20.d0*u(1) &
	+15.d0*u(2)-6.d0*u(3)+u(4))/h6
	d6f(2)=(u(-1)-6.d0*u(0)+15.d0*u(1)-20.d0*u(2) &
	+15.d0*u(3)-6.d0*u(4)+u(5))/h6

	r(1)=(3.d0*sig*u(3)+(42.d0+2.d0*sig)*u(2)+(-48.d0+12.d0*sig)*u(1) &
	+(6.d0-18.d0*sig)*u(0)+sig*u(-1))/12.d0/h-(1.d0-sig)*h*s(1) &
	& +(2.d0+sig)*esi*h5*d6f(1)+(1.d0+sig)*esi*h5*d6f(2)
	r(2)=(-3.d0*sig*u(0)+(-42.d0-2.d0*sig)*u(1)+(48.d0-12.d0*sig)*u(2) &
	+(-6.d0+18.d0*sig)*u(3)-sig*u(4))/12.d0/h+(1.d0-sig)*h*s(2)

	f(1)=((2.d0+sig)*r(1)-(1.d0+sig)*r(2))/(3.d0+2.d0*sig)
end subroutine

!---------------------------------------------------------------------

subroutine OCFD_DF_SB2_UCC45_M(u,f,s,n,h)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h,sig,esi,h2,h4,h5,h6
    integer n
	real(kind=OCFD_REAL_KIND)::r(n-1:n),d6f(n-1:n)
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),f(1-LAP:n+LAP)
	real(kind=OCFD_REAL_KIND)::s(1-LAP:n+LAP)

	sig=0.54395d0
!	esi=-0.01385977393639451d0
	esi=-0.01820010801241671d0
	h2=h*h
	h4=h2*h2
	h5=h4*h
	h6=h4*h2

	d6f(n)=(u(n-3)-6.d0*u(n-2)+15.d0*u(n-1)-20.d0*u(n) &
	+15.d0*u(n+1)-6.d0*u(n+2)+u(n+3))/h6
	d6f(n-1)=(u(n-4)-6.d0*u(n-3)+15.d0*u(n-2)-20.d0*u(n-1) &
	+15.d0*u(n)-6.d0*u(n+1)+u(n+2))/h6

	r(n-1)=(3.d0*sig*u(n+1)+(42.d0+2.d0*sig)*u(n)+(-48.d0+12.d0*sig)*u(n-1) &
	+(6.d0-18.d0*sig)*u(n-2)+sig*u(n-3))/12.d0/h-(1.d0-sig)*h*s(n-1) 
	r(n)=(-3.d0*sig*u(n-2)+(-42.d0-2.d0*sig)*u(n-1)+(48.d0-12.d0*sig)*u(n) &
	+(-6.d0+18.d0*sig)*u(n+1)-sig*u(n+2))/12.d0/h+(1.d0-sig)*h*s(n) &
	& -(2.d0+sig)*esi*h5*d6f(n)-(1.d0+sig)*esi*h5*d6f(n-1)

	f(n)=((2.d0+sig)*r(n)-(1.d0+sig)*r(n-1))/(3.d0+2.d0*sig)
end subroutine

!!!================================================================!!!






!!!=======================D2F_BOUND_PADE4=========================!!!
!!boundary
subroutine OCFD_D2F_BOUND_PADE4(u,s,n,h,flag)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h,h2,pi,pi2
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),s(1-LAP:n+LAP)
	integer n,flag

	pi=4.d0*datan(1.d0)
	pi2=pi*pi
	h2=h*h

	if (flag.eq.1) then
!		s(1)=-4.d0*pi2*dsin(2.d0*pi*0.d0)-4.d0*pi2*dcos(2.d0*pi*0.d0)
!		s(1)=(4.d0*(0.d0-5)**2-2.d0)*exp(-(0.d0-5.d0)**2)
		s(1)=0.d0
	elseif (flag.eq.2) then
!		s(n)=-4.d0*pi2*dsin(2.d0*pi*1.d0)-4.d0*pi2*dcos(2.d0*pi*1.d0)
!		s(n)=(4.d0*(10.d0-5)**2-2.d0)*exp(-(10.d0-5.d0)**2)
		s(n)=0.d0
	endif

end subroutine


!! Subboundary
subroutine OCFD_D2F_SB_PADE4(u,s,n,h,flag)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h,h2
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),s(1-LAP:n+LAP)
	integer n,flag

	h2=h*h

	if (flag.eq.1) then
		s(n)=((1.d0/144.d0)*u(n-3)-(1.d0/8.d0)*u(n-2)+(23.d0/16.d0)*u(n-1)-(95.d0/36.d0)*u(n) &
		+(23.d0/16.d0)*u(n+1)-(1.d0/8.d0)*u(n+2)+(1.d0/144.d0)*u(n+3))/h2
	elseif (flag.eq.2) then
		s(1)=((1.d0/144.d0)*u(-2)-(1.d0/8.d0)*u(-1)+(23.d0/16.d0)*u(0)-(95.d0/36.d0)*u(1) &
		+(23.d0/16.d0)*u(2)-(1.d0/8.d0)*u(3)+(1.d0/144.d0)*u(4))/h2
	elseif (flag.eq.0) then
		s(1)=((1.d0/144.d0)*u(-2)-(1.d0/8.d0)*u(-1)+(23.d0/16.d0)*u(0)-(95.d0/36.d0)*u(1) &
		+(23.d0/16.d0)*u(2)-(1.d0/8.d0)*u(3)+(1.d0/144.d0)*u(4))/h2
		s(n)=((1.d0/144.d0)*u(n-3)-(1.d0/8.d0)*u(n-2)+(23.d0/16.d0)*u(n-1)-(95.d0/36.d0)*u(n) &
		+(23.d0/16.d0)*u(n+1)-(1.d0/8.d0)*u(n+2)+(1.d0/144.d0)*u(n+3))/h2
	endif
end subroutine

!!!============================================================!!!


!!!=========================D6F_BOUND_CENTER=========================!!!
!boundary
subroutine OCFD_D6F_BOUND_CENTER(u,d6f,n,h,flag)
	include 'openNS3d.h'		
	real(kind=OCFD_REAL_KIND)::h,h2,h4,h6
	real(kind=OCFD_REAL_KIND)::u(1-LAP:n+LAP),d6f(1-LAP:n+LAP)
	integer n,flag

	h2=h*h
	h4=h2*h2
	h6=h4*h2

	if (flag==1) then
		d6f(1)=0.d0
		d6f(2)=0.d0
		d6f(3)=0.d0
	elseif (flag==2) then
		d6f(n)=0.d0
		d6f(n-1)=0.d0
		d6f(n-2)=0.d0
	endif

end subroutine
