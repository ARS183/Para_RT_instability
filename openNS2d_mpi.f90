subroutine part2d()	

	include 'openNS3d.h'

	integer :: i,j,k,ka,npx1,npy1,npx2,npy2,my_mod1,i_global,j_global
!	real(kind=OCFD_REAL_KIND) :: hx,hy
!	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: xx,yy
    real(kind=OCFD_REAL_KIND) :: pi !!

    pi=datan(1.d0)*4.d0 !!
	hx=slx/dble(nx_global-1)!!
	hy=sly/dble(ny_global-1)!!

	if (np_size .ne. npx0*npy0) then	
		if (my_id .eq. 0) then
			write(*,*)'Number of Total Procs is ',np_size
		  		write(*,"(' Number of Total procs is not equal to npx0*npy0=',I3,' !!!')")npx0*npy0
		endif
		stop
	endif


	npx=mod(my_id,npx0)
	npy=my_id/npx0

	CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,   npy, npx,MPI_COMM_X,ierr)		! in a x-line 
	CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,   npx, npy,MPI_COMM_Y,ierr)		! in a y-line

	nx=nx_global/npx0
	ny=ny_global/npy0
	if(npx .lt. mod(nx_global,npx0)) nx=nx+1
	if(npy .lt. mod(ny_global,npy0)) ny=ny+1

	do k=0,npx0-1
		ka=min(k,mod(nx_global,npx0))
		i_offset(k)=int(nx_global/npx0)*k+ka+1
		i_nn(k)=nx_global/npx0
		 	if(k .lt. mod(nx_global,npx0)) i_nn(k)=i_nn(k)+1
	enddo

	do k=0,npy0-1
		ka=min(k,mod(ny_global,npy0))
		j_offset(k)=int(ny_global/npy0)*k+ka+1
		j_nn(k)=ny_global/npy0
		 	if(k .lt. mod(ny_global,npy0)) j_nn(k)=j_nn(k)+1
	enddo
	
	

	npx1=my_mod1(npx-1,npx0)
	npx2=my_mod1(npx+1,npx0)
	ID_XM1=npy*npx0+npx1		! -1 proc in x-direction
	ID_XP1=npy*npx0+npx2		! +1 proc in x-direction
	if(Iperiodic_X .eq.0 .and. npx .eq. 0) ID_XM1=MPI_PROC_NULL		 ! if not periodic, 0 node donot send mesg to npx0-1 node
	if(Iperiodic_X .eq.0 .and. npx .eq. npx0-1) ID_XP1=MPI_PROC_NULL
		  
	npy1=my_mod1(npy-1,npy0)
	npy2=my_mod1(npy+1,npy0)
	ID_YM1=npy1*npx0+npx
	ID_YP1=npy2*npx0+npx
	if(Iperiodic_Y .eq.0 .and. npy .eq. 0) ID_YM1=MPI_PROC_NULL		 ! if not periodic, 0 node donot send mesg to npy0-1 node
	if(Iperiodic_Y .eq.0 .and. npy .eq. npy0-1) ID_YP1=MPI_PROC_NULL


end subroutine




	function my_mod1(i,n)
	include 'openNS3d.h'
	integer :: my_mod1,i,n
		 	if(i.lt.0) then
		 	my_mod1=i+n
		 	else if (i.gt.n-1) then
		 	my_mod1=i-n
		 	else
		 	my_mod1=i
		 	endif
	end


	subroutine check_xy2d(f)
		include 'openNS3d.h'
		real(kind=OCFD_REAL_KIND) :: f(1-LAP:nx+LAP,1-LAP:ny+LAP)
        !        integer ierr                   !!no2-2
		call check_x2d(f)
        !        call MPI_BARRIER(MPI_COMM_WORLD,ierr)   !!no2-2
		call check_y2d(f)
		return
	end subroutine


	subroutine check_x2d(f)
		include 'openNS3d.h'
!send and recv mesg, to check array in x direction.
		
		integer i,j,k1,npx1,npx2,my_mod1
		real(kind=OCFD_REAL_KIND):: f(1-LAP:nx+LAP,1-LAP:ny+LAP),   &
		    tmp_send1(LAP*ny),tmp_send2(LAP*ny), tmp_recv1(LAP*ny),tmp_recv2(LAP*ny)
	  
		
		do j=1,ny
		do i=1,LAP
		k1=(j-1)*LAP+i
		tmp_send1(k1)=f(i,j)
		tmp_send2(k1)=f(nx-LAP+i,j)
		enddo
		enddo
!!
!               npx=mod(my_id,npx0)
!               npy=my_id/npx0      
               if (npx .eq. 0) then
               do j=1,ny
               do i=1,LAP
               k1=(j-1)*LAP+i
               tmp_send1(k1)=f(i+1,j)
               enddo
               enddo
               endif

               if (npx .eq. npx0-1) then
               do j=1,ny
               do i=1,LAP
               k1=(j-1)*LAP+i
               tmp_send2(k1)=f(nx-LAP+i-1,j)
               enddo
               enddo
               endif
!!

		call MPI_Sendrecv(tmp_send1,LAP*ny,  OCFD_DATA_TYPE, ID_XM1, 9000, &
		    tmp_recv2, LAP*ny,  OCFD_DATA_TYPE,ID_XP1, 9000,MPI_COMM_WORLD,Status,ierr)
		call MPI_Sendrecv(tmp_send2,LAP*ny,  OCFD_DATA_TYPE,ID_XP1,8000,   &
		    tmp_recv1,LAP*ny, OCFD_DATA_TYPE,ID_XM1,  8000,MPI_COMM_WORLD,Status,ierr)

		if(ID_XM1 .ne. MPI_PROC_NULL) then
		 do j=1,ny
		 do i=1,LAP
		   k1=(j-1)*LAP+i
		   f(i-LAP,j)=tmp_recv1(k1)
		 enddo
		 enddo
		endif

		if(ID_XP1 .ne. MPI_PROC_NULL) then
		 do j=1,ny
		 do i=1,LAP
		   k1=(j-1)*LAP+i
		   f(nx+i,j)=tmp_recv2(k1)
		 enddo
		 enddo
		endif

	end subroutine



	subroutine check_y2d(f)
		include 'openNS3d.h'
		
		integer i,j,k1,npy1,npy2,my_mod1
		real(kind=OCFD_REAL_KIND):: f(1-LAP:nx+LAP,1-LAP:ny+LAP),  &
		    tmp_send1(LAP*(nx+2*LAP)),tmp_send2(LAP*(nx+2*LAP)), &
                      tmp_recv1(LAP*(nx+2*LAP)),tmp_recv2(LAP*(nx+2*LAP))    !!no2

		do j=1,LAP
		do i=1-LAP,nx+LAP      !!no2
		k1=(j-1)*(nx+2*LAP)+i+LAP
		tmp_send1(k1)=f(i,j)
		tmp_send2(k1)=f(i,ny+j-LAP)
		enddo
		enddo
!!
!                npx=mod(my_id,npx0)
!               npy=my_id/npx0
                if (npy .eq. 0) then
		do j=1,LAP
		do i=1-LAP,nx+LAP   !!no2
		k1=(j-1)*(nx+2*LAP)+i+LAP
		tmp_send1(k1)=f(i,j+1)
!		tmp_send2(k1)=f(i,ny+j-LAP)
		enddo
                enddo
                endif

                if (npy .eq. npy0-1) then
		do j=1,LAP
		do i=1-LAP,nx+LAP  !!no2
		k1=(j-1)*(nx+2*LAP)+i+LAP
!		tmp_send1(k1)=f(i,j)
		tmp_send2(k1)=f(i,ny-LAP+j-1)
		enddo
                enddo
                endif
!!
		call MPI_Sendrecv(tmp_send1,LAP*(nx+2*LAP),  OCFD_DATA_TYPE, ID_YM1, 9000, &
		    tmp_recv2, LAP*(nx+2*LAP),  OCFD_DATA_TYPE,ID_YP1, 9000,MPI_COMM_WORLD,Status,ierr)   !!no2
		call MPI_Sendrecv(tmp_send2,LAP*(nx+2*LAP),  OCFD_DATA_TYPE,ID_YP1,8000,   &
		    tmp_recv1,LAP*(nx+2*LAP), OCFD_DATA_TYPE,ID_YM1,  8000,MPI_COMM_WORLD,Status,ierr)    !!no2


		if(ID_YM1 .ne. MPI_PROC_NULL) then
		 do j=1,LAP
		 do i=1-LAP,nx+LAP    !!no2
		   k1=(j-1)*(nx+2*LAP)+i+LAP
		   f(i,j-LAP)=tmp_recv1(k1)
		 enddo
		 enddo
		endif

		if(ID_YP1 .ne. MPI_PROC_NULL) then
		 do j=1,LAP
		 do i=1-LAP,nx+LAP    !!no2
		   k1=(j-1)*(nx+2*LAP)+i+LAP
		   f(i,ny+j)=tmp_recv2(k1)
		 enddo
		 enddo
		endif

	end subroutine



