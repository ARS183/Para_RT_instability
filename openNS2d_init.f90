
!!!==========================================================================================
!!!Init Taylor Vortex
!!!==========================================================================================

subroutine init_TV(rho,u,v,p,xx,yy) !!
	
	include 'openNS3d.h'
	integer :: i,j,i_global,j_global

	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: u,v,p,xx,yy,rho
!    real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: uacc,vacc,pacc !!
	real(kind=OCFD_REAL_KIND) :: parat,spacc,rou,sigma   !!!
	real(kind=OCFD_REAL_KIND) :: tempexp,xbeg,ybeg
	integer :: paran  !!!
	
	real(kind=OCFD_REAL_KIND),parameter::xo=4.d0,yo=0.d0
	real(kind=OCFD_REAL_KIND),parameter::M0=0.4d0,Mv=0.2d0,gama=1.4d0,p0=1.d0/gama/M0**2.d0,rho0=1.d0,u0=1.d0
	real(kind=OCFD_REAL_KIND)::Ms=M0,ps=p0

!

	xbeg=0.d0
	ybeg=-4.d0

	do i=1,nx
	do j=1,ny
	
		i_global=i_offset(npx)+i-1
		j_global=j_offset(npy)+j-1
		xx(i,j)=xbeg+dble(i_global-1)*hx
		yy(i,j)=ybeg+dble(j_global-1)*hy

		tempexp=dexp((1.d0-(xx(i,j)-xo)**2.d0-(yy(i,j)-yo)**2.d0)/2.d0)

		p(i,j)=1.d0/gama/M0**2.d0*(1.d0-(gama-1.d0)/2.d0*Mv**2.d0*tempexp**2.d0)**(gama/(gama-1.d0))
		rho(i,j)=(1.d0-(gama-1.d0)/2.d0*Mv**2.d0*tempexp**2.d0)**(1.d0/(gama-1.d0))
		u(i,j)=-(yy(i,j)-yo)*Mv/M0*tempexp+u0
		v(i,j)=(xx(i,j)-xo)*Mv/M0*tempexp

	enddo
	enddo
	
!	if (npx==0) then
!		p(1,:)=p0
!		rho(1,:)=rho0
!		u(1,:)=u0
!		v(1,:)=0.d0
!	endif

!	if (npx .eq. 0) then
!		u(1,:)=0.d0
!		v(1,:)=0.d0
!	endif

!!
	end subroutine
