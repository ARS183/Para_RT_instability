	subroutine GS_solver(u,f,istep,hx,hy,NUM_METHOD)
	
	include 'openNS2d.h'

	integer :: i,j,k,istep,nx1,nx2,ny1,ny2,npx2,npy2,NUM_METHOD
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: f,u,u_old,s
	real(kind=OCFD_REAL_KIND),dimension(1:nx,1:ny) :: ap,as,an,ae,aw,ase,asw,ane,anw
	real(kind=OCFD_REAL_KIND) :: error,error1,hx,hy
	real(kind=OCFD_REAL_KIND),dimension(ny*LAP) :: sendx,recvx,sendx1
	real(kind=OCFD_REAL_KIND),dimension(nx*LAP) :: sendy,recvy,sendy1
        real(kind=OCFD_REAL_KIND) :: pi

        pi=datan(1.d0)*4.d0 !!
	u_old=0.d0

	nx1=1
	nx2=nx
	ny1=1
	ny2=ny
        if (Iperiodic_X .ne. 1) then !!
	if (npx .eq. 0) then
		nx1=2
	endif
	if (npx .eq. npx0-1) then
		nx2=nx-1
	endif
        endif !!
        if (Iperiodic_Y .ne. 1) then !!
	if (npy .eq. 0) then
		ny1=2
	endif
	if (npy .eq. npy0-1) then
		ny2=ny-1
	endif
        endif !!

	if (NUM_METHOD .eq. 1) then
	call coeff_poisson_2nd(ap,as,an,ae,aw,ase,asw,ane,anw,f,s,hx,hy)
	elseif (NUM_METHOD .eq. 2) then
	call coeff_poisson_4th(ap,as,an,ae,aw,ase,asw,ane,anw,f,s,hx,hy)
	elseif (NUM_METHOD .eq. 3) then
	call coeff_poisson_6th(ap,as,an,ae,aw,ase,asw,ane,anw,f,s,hx,hy)
	endif
	error=1.d0
	istep=0


	do while(error .gt. eps)

		istep=istep+1
			
		u_old=u		
	
		npx2=(npx0+1)/2
		npy2=(npy0+1)/2
		
		if (npx.lt.npx2 .and. npy.lt.npy2) then

		do i=nx1,nx2
		do j=ny1,ny2
		u(i,j)=(s(i,j)-ae(i,j)*u(i+1,j)-aw(i,j)*u(i-1,j)-an(i,j)*u(i,j+1)-as(i,j)*u(i,j-1)-ase(i,j)*u(i+1,j-1) &
		-asw(i,j)*u(i-1,j-1)-ane(i,j)*u(i+1,j+1)-anw(i,j)*u(i-1,j+1))/ap(i,j)*rf+(1.d0-rf)*u(i,j)
		enddo
		enddo

	
		elseif (npx.ge.npx2 .and. npy.lt.npy2) then

		do i=nx2,nx1,-1
		do j=ny1,ny2
		u(i,j)=(s(i,j)-ae(i,j)*u(i+1,j)-aw(i,j)*u(i-1,j)-an(i,j)*u(i,j+1)-as(i,j)*u(i,j-1)-ase(i,j)*u(i+1,j-1) &
		-asw(i,j)*u(i-1,j-1)-ane(i,j)*u(i+1,j+1)-anw(i,j)*u(i-1,j+1))/ap(i,j)*rf+(1.d0-rf)*u(i,j)
		enddo
		enddo

		elseif (npx.lt.npx2 .and. npy.ge.npy2) then

		do i=nx1,nx2
		do j=ny2,ny1,-1
		u(i,j)=(s(i,j)-ae(i,j)*u(i+1,j)-aw(i,j)*u(i-1,j)-an(i,j)*u(i,j+1)-as(i,j)*u(i,j-1)-ase(i,j)*u(i+1,j-1) &
		-asw(i,j)*u(i-1,j-1)-ane(i,j)*u(i+1,j+1)-anw(i,j)*u(i-1,j+1))/ap(i,j)*rf+(1.d0-rf)*u(i,j)
		enddo
		enddo

		else

		do i=nx2,nx1,-1
		do j=ny2,ny1,-1
		u(i,j)=(s(i,j)-ae(i,j)*u(i+1,j)-aw(i,j)*u(i-1,j)-an(i,j)*u(i,j+1)-as(i,j)*u(i,j-1)-ase(i,j)*u(i+1,j-1) &
		-asw(i,j)*u(i-1,j-1)-ane(i,j)*u(i+1,j+1)-anw(i,j)*u(i-1,j+1))/ap(i,j)*rf+(1.d0-rf)*u(i,j)
		enddo
		enddo

		endif

!===only for pressure poisson equation
                
		if (NUM_METHOD .eq. 1) then
		call poisson_2nd_BC(u,f,hx,hy)
		elseif (NUM_METHOD .eq. 2) then
		call poisson_4th_BC(u,f,hx,hy)
		elseif (NUM_METHOD .eq. 3) then
		call poisson_6th_BC(u,f,hx,hy)
		endif
		call poisson_average(u)

		error1=maxval(dabs(u(1:nx,1:ny)-u_old(1:nx,1:ny)))
	
		call MPI_ALLREDUCE(error1,error,1,OCFD_DATA_TYPE,MPI_MAX,MPI_COMM_WORLD,ierr)

!		if (my_id .eq. 0) then
!			write(*,*)Istep,error
!		endif

		call check_xy2d(u)

	end do ! do while


!	if (my_id .eq. 0) then
!		write(*,*)'================================='
!		write(*,*)'Congratulation!!! It''s finished.'
!		write(*,"(' Poisson Eqn Istep=',I5,' Max error=',E16.6)")Istep,error
!		write(*,*)'================================='
!	endif

	end subroutine



	subroutine coeff_poisson_2nd(ap,as,an,ae,aw,ase,asw,ane,anw,f,s,hx,hy)
	include 'openNS2d.h'

	integer :: i,j
	real(kind=OCFD_REAL_KIND) :: hx,hy,hx2,hy2
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP)::f,s
	real(kind=OCFD_REAL_KIND),dimension(1:nx,1:ny) :: ap,as,an,ae,aw,ase,asw,ane,anw

	hx2=hx*hx
	hy2=hy*hy

	do 100 j=1,ny
	do 100 i=1,nx		
		ap(i,j)=-2.d0*(1.d0/hx2+1.d0/hy2)
		ae(i,j)=1.d0/hx2
		aw(i,j)=1.d0/hx2
		an(i,j)=1.d0/hy2
		as(i,j)=1.d0/hy2
		ase(i,j)=0.d0
		asw(i,j)=0.d0
		ane(i,j)=0.d0
		anw(i,j)=0.d0
		s(i,j)=f(i,j)
100	continue

	end subroutine
	
	subroutine coeff_poisson_4th(ap,as,an,ae,aw,ase,asw,ane,anw,f,s,hx,hy)
	include 'openNS2d.h'

	integer :: i,j
	real(kind=OCFD_REAL_KIND) :: hx,hy,hx2,hy2
        real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP)::f,s
	real(kind=OCFD_REAL_KIND),dimension(1:nx,1:ny) :: ap,as,an,ae,aw,ase,asw,ane,anw

	hx2=hx*hx
	hy2=hy*hy

	do 100 j=1,ny
	do 100 i=1,nx		
		ap(i,j)=-20.d0*(1.d0/hx2+1.d0/hy2)
		ae(i,j)=2.d0*(5.d0/hy2-1.d0/hx2)
		aw(i,j)=2.d0*(5.d0/hy2-1.d0/hx2)
		an(i,j)=2.d0*(5.d0/hx2-1.d0/hy2)
		as(i,j)=2.d0*(5.d0/hx2-1.d0/hy2)
		ase(i,j)=(1.d0/hy2+1.d0/hx2)
		asw(i,j)=(1.d0/hy2+1.d0/hx2)
		ane(i,j)=(1.d0/hy2+1.d0/hx2)
		anw(i,j)=(1.d0/hy2+1.d0/hx2)
		s(i,j)=(8.d0*f(i,j)+f(i-1,j)+f(i,j+1)+f(i+1,j)+f(i,j-1))
100	continue

	end subroutine	


	subroutine coeff_poisson_6th(ap,as,an,ae,aw,ase,asw,ane,anw,f,s,hx,hy)
	include 'openNS2d.h'

	integer :: i,j,flag
	real(kind=OCFD_REAL_KIND) :: hx,hy,hx2,hy2
        real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP)::f,s,fxx,fyy
	real(kind=OCFD_REAL_KIND),dimension(1:nx,1:ny) :: ap,as,an,ae,aw,ase,asw,ane,anw
	
	flag=1

	hx2=hx*hx
	hy2=hy*hy
	

	
	do i=1,nx
	call OCFD_D2F_SCD(f(i,:),fyy(i,:),ny,hy,4,flag)
	enddo

	do j=1,ny
	call OCFD_D2F_SCD(f(:,j),fxx(:,j),nx,hx,4,flag)
	enddo
    	
    	
	
	do 100 j=1,ny
	do 100 i=1,nx		
		ap(i,j)=-10.d0*(1.d0/hx2+1.d0/hy2)
		ae(i,j)=4.d0/hx2
		aw(i,j)=4.d0/hx2
		an(i,j)=4.d0/hy2
		as(i,j)=4.d0/hy2
		ase(i,j)=1.d0/hy2
		asw(i,j)=1.d0/hy2
		ane(i,j)=1.d0/hy2
		anw(i,j)=1.d0/hy2
		s(i,j)=1d0/15d0*(82.d0*f(i,j)+f(i-1,j)+f(i,j+1)+f(i+1,j)+f(i,j-1) &
			+f(i-1,j-1)+f(i-1,j+1)+f(i+1,j+1)+f(i+1,j-1)) &
			+0.3d0*(hx2*fxx(i,j)+hy2*fyy(i,j))
100	continue

	end subroutine	


	subroutine poisson_2nd_BC(p,f,hx,hy)
	include 'openNS2d.h'
	integer :: i,j
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: p,f
	real(kind=OCFD_REAL_KIND) :: hx,hy,hx2,hy2
!	write(*,"('==',I4,4E14.4)")my_id,maxval(p),minval(p)

	if (Iperiodic_X .ne. 1) then
	if (npx.eq.0) then
        i=1
        do j=1,ny
        p(i,j)=p(i+1,j)
        enddo
	endif

    	if (npx.eq.npx0-1) then
        i=nx
        do j=1,ny
        p(i,j)=p(i-1,j)
        enddo
	endif

	endif

	if (Iperiodic_Y .ne. 1) then
    	if (npy.eq.0) then
        j=1
        do i=1,nx
        p(i,j)=p(i,j+1)
        enddo
	endif

    	if (npy.eq.npy0-1) then
        j=ny
        do i=1,nx
        p(i,j)=p(i,j-1)
        enddo
	endif
	
	endif

	if (Iperiodic_X .ne. 1 .and. Iperiodic_Y .ne. 1) then

	if (npx.eq.0 .and. npy.eq.0) then
        i=1
        j=1
        p(i,j)=(p(i,j+1)+p(i+1,j))/2.d0
    	endif

    	if (npx.eq.0 .and. npy.eq.npy0-1) then
        i=1
        j=ny
        p(i,j)=(p(i,j-1)+p(i+1,j))/2.d0
    	endif

    	if (npx.eq.npx0-1 .and. npy.eq.0) then
        i=nx
        j=1
        p(i,j)=(p(i,j+1)+p(i-1,j))/2.d0
    	endif
	
    	if (npx.eq.npx0-1 .and. npy.eq.npy0-1) then
        i=nx
        j=ny
        p(i,j)=(p(i,j-1)+p(i-1,j))/2.d0
    	endif

	endif

!	write(*,"('oo',I4,4E14.4)")my_id,maxval(p),minval(p)
	end subroutine
	
	
	subroutine poisson_4th_BC(p,f,hx,hy)
	include 'openNS2d.h'
	integer :: i,j
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: p,f
	real(kind=OCFD_REAL_KIND) :: hx,hy,hx2,hy2
!	write(*,"('==',I4,4E14.4)")my_id,maxval(p),minval(p)

	hx2=hx*hx
	hy2=hy*hy
	

	if (Iperiodic_X .ne. 1) then
	if (npx.eq.0) then
        i=1
        do j=1,ny
        p(i,j)=(p(i+1,j)/hx2+0.5d0*(p(i,j+1)+p(i,j-1))/hy2-0.5d0*f(i,j))/(1.d0/hx2+1.d0/hy2)
        enddo
	endif

    	if (npx.eq.npx0-1) then
        i=nx
        do j=1,ny
        p(i,j)=(p(i-1,j)/hx2+0.5d0*(p(i,j+1)+p(i,j-1))/hy2-0.5d0*f(i,j))/(1.d0/hx2+1.d0/hy2)
        enddo
	endif

	endif

	if (Iperiodic_Y .ne. 1) then
    	if (npy.eq.0) then
        j=1
        do i=1,nx
         p(i,j)=(p(i,j+1)/hy2+0.5d0*(p(i+1,j)+p(i-1,j))/hx2-0.5d0*f(i,j))/(1.d0/hx2+1.d0/hy2)
        enddo
	endif

    	if (npy.eq.npy0-1) then
        j=ny
        do i=1,nx
        p(i,j)=(p(i,j-1)/hy2+0.5d0*(p(i+1,j)+p(i-1,j))/hx2-0.5d0*f(i,j))/(1.d0/hx2+1.d0/hy2)
        enddo
	endif
	
	endif

	if (Iperiodic_X .ne. 1 .and. Iperiodic_Y .ne. 1) then

	if (npx.eq.0 .and. npy.eq.0) then
        i=1
        j=1
         p(i,j)=(p(i+1,j)/hx2+p(i,j+1)/hy2-0.5d0*f(i,j))/(1.d0/hx2+1.d0/hy2)
    	endif

    	if (npx.eq.0 .and. npy.eq.npy0-1) then
        i=1
        j=ny
        p(i,j)=(p(i+1,j)/hx2+p(i,j-1)/hy2-0.5d0*f(i,j))/(1.d0/hx2+1.d0/hy2)
    	endif

    	if (npx.eq.npx0-1 .and. npy.eq.0) then
        i=nx
        j=1
        p(i,j)=(p(i-1,j)/hx2+p(i,j+1)/hy2-0.5d0*f(i,j))/(1.d0/hx2+1.d0/hy2)
    	endif
	
    	if (npx.eq.npx0-1 .and. npy.eq.npy0-1) then
        i=nx
        j=ny
        p(i,j)=(p(i-1,j)/hx2+p(i,j-1)/hy2-0.5d0*f(i,j))/(1.d0/hx2+1.d0/hy2)
    	endif

	endif

!	write(*,"('oo',I4,4E14.4)")my_id,maxval(p),minval(p)
	end subroutine




	subroutine poisson_6th_BC(p,f,hx,hy)
	include 'openNS2d.h'
	integer :: i,j
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: p,f
	real(kind=OCFD_REAL_KIND) :: hx,hy,hx2,hy2
!	write(*,"('==',I4,4E14.4)")my_id,maxval(p),minval(p)

	if (Iperiodic_X .ne. 1) then
	if (npx.eq.0) then
        i=1
        do j=1,ny
        p(i,j)=p(i+1,j)
        enddo
	endif

    	if (npx.eq.npx0-1) then
        i=nx
        do j=1,ny
        p(i,j)=p(i-1,j)
        enddo
	endif

	endif

	if (Iperiodic_Y .ne. 1) then
    	if (npy.eq.0) then
        j=1
        do i=1,nx
        p(i,j)=p(i,j+1)
        enddo
	endif

    	if (npy.eq.npy0-1) then
        j=ny
        do i=1,nx
        p(i,j)=p(i,j-1)
        enddo
	endif
	
	endif

	if (Iperiodic_X .ne. 1 .and. Iperiodic_Y .ne. 1) then

	if (npx.eq.0 .and. npy.eq.0) then
        i=1
        j=1
        p(i,j)=(p(i,j+1)+p(i+1,j))/2.d0
    	endif

    	if (npx.eq.0 .and. npy.eq.npy0-1) then
        i=1
        j=ny
        p(i,j)=(p(i,j-1)+p(i+1,j))/2.d0
    	endif

    	if (npx.eq.npx0-1 .and. npy.eq.0) then
        i=nx
        j=1
        p(i,j)=(p(i,j+1)+p(i-1,j))/2.d0
    	endif
	
    	if (npx.eq.npx0-1 .and. npy.eq.npy0-1) then
        i=nx
        j=ny
        p(i,j)=(p(i,j-1)+p(i-1,j))/2.d0
    	endif

	endif

!	write(*,"('oo',I4,4E14.4)")my_id,maxval(p),minval(p)
	end subroutine


	subroutine poisson_average(p)
	include 'openNS2d.h'
	integer :: i,j
	real(kind=OCFD_REAL_KIND) :: psum,psum_total
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: p

	psum=0.d0
	do i=1,nx
	do j=1,ny
	psum=psum+p(i,j)
	enddo
	enddo
	call MPI_ALLREDUCE(psum,psum_total,1,OCFD_DATA_TYPE,MPI_SUM,MPI_COMM_WORLD,ierr)
!	psum_total=psum/dble(nx_global)/dble(ny_global)
	psum_total=psum_total/dble(nx_global)/dble(ny_global)
	
	p=p-psum_total	

!	do j=1,ny
!   do i=1,nx
!	p(i,j)=p(i,j)-psum_total
!   err=max(err,dabs(p00(i,j)-p(i,j)))
!   enddo
!	enddo

!	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!	call check_xy2d(p)


	end subroutine


!	subroutine coeff_2CDS(ap,as,an,ae,aw,hx,hy)
!	include 'openNS2d.h'
!
!	integer :: i,j
!	real(kind=OCFD_REAL_KIND) :: y,hx,hy
!	real(kind=OCFD_REAL_KIND),dimension(1:nx,1:ny) :: ap,as,an,ae,aw
!
!	do 100 j=1,ny
!		y=(j_offset(npy)+j-2)*hy
!	do 101 i=1,nx		
!		ap(i,j)=2.d0*eps*(1.d0/hx**2.d0+1.d0/hy**2.d0)
!		ae(i,j)=-eps/hx**2.d0
!		aw(i,j)=-eps/hx**2.d0
!		an(i,j)=-eps/hy**2.d0+1.d0/(1.d0+y)/2.d0/hy
!		as(i,j)=-eps/hy**2.d0-1.d0/(1.d0+y)/2.d0/hy
!101	continue
!100	continue
!
!	end subroutine

	
