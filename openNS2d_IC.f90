!===BCs for different problems

	subroutine CASE_LDC(u,v,hx,hy)

	include 'openNS2d.h'

	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: u,v
	real(kind=OCFD_REAL_KIND) :: hx,hy

	if (npx .eq. 0) then
		u(1,:)=0.d0
		v(1,:)=0.d0
	endif

	
	if (npx .eq. npx0-1) then
		u(nx,:)=0.d0
		v(nx,:)=0.d0
	endif

	if (npy .eq. 0) then
		u(:,1)=0.d0
		v(:,1)=0.d0
	endif

	if (npy .eq. npy0-1) then
		u(:,ny)=1.d0
		v(:,ny)=0.d0
	endif

	end subroutine
