	subroutine write_data(rho,u,v,p,xx,yy,filename)
	!subroutine write_data(rho,u,v,p,xx,yy)

	include 'openNS3d.h'

	integer :: i,j,k,i_global,j_global,ii,jj
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP,1-LAP:ny+LAP) :: rho,u,v,p,xx,yy
	real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:) :: temp_rho,temp_u,temp_v
	real(kind=OCFD_REAL_KIND),allocatable,dimension(:,:) :: temp_p,temp_xx,temp_yy
	real(kind=OCFD_REAL_KIND),dimension(1:nx_global,1:ny_global) :: rho_global,u_global,v_global
	real(kind=OCFD_REAL_KIND),dimension(1:nx_global,1:ny_global) :: p_global,xx_global,yy_global

		
	if (my_id .eq. 0 ) then
		do 100 j=0,npy0-1
		do 100 i=0,npx0-1		
			k=j*npx0+i
			
			if (k .eq. 0) then
				do jj=1,ny
				do ii=1,nx
					i_global=i_offset(i)-1+ii
					j_global=j_offset(j)-1+jj
					rho_global(i_global,j_global)=rho(ii,jj)
					u_global(i_global,j_global)=u(ii,jj)
					v_global(i_global,j_global)=v(ii,jj)
					p_global(i_global,j_global)=p(ii,jj)					
					xx_global(i_global,j_global)=xx(ii,jj)
					yy_global(i_global,j_global)=yy(ii,jj)
				enddo
				enddo
			else
				allocate(temp_u(i_nn(i),j_nn(j)),temp_v(i_nn(i),j_nn(j)),temp_p(i_nn(i),j_nn(j)), &
					temp_xx(i_nn(i),j_nn(j)),temp_yy(i_nn(i),j_nn(j)),temp_rho(i_nn(i),j_nn(j)))
				
				call MPI_RECV(temp_rho,i_nn(i)*j_nn(j),OCFD_DATA_TYPE,  &
				k,6,MPI_COMM_WORLD,status,ierr) 
				call MPI_RECV(temp_u,i_nn(i)*j_nn(j),OCFD_DATA_TYPE,  &
				k,1,MPI_COMM_WORLD,status,ierr) 
				call MPI_RECV(temp_v,i_nn(i)*j_nn(j),OCFD_DATA_TYPE,  &
				k,2,MPI_COMM_WORLD,status,ierr) 
				call MPI_RECV(temp_p,i_nn(i)*j_nn(j),OCFD_DATA_TYPE,  &
				k,3,MPI_COMM_WORLD,status,ierr)
				call MPI_RECV(temp_xx,i_nn(i)*j_nn(j),OCFD_DATA_TYPE,  &
				k,4,MPI_COMM_WORLD,status,ierr) 
				call MPI_RECV(temp_yy,i_nn(i)*j_nn(j),OCFD_DATA_TYPE,  &
				k,5,MPI_COMM_WORLD,status,ierr)

				do jj=1,j_nn(j)
				do ii=1,i_nn(i)
					i_global=i_offset(i)+ii-1
					j_global=j_offset(j)+jj-1
					rho_global(i_global,j_global)=temp_rho(ii,jj)
					u_global(i_global,j_global)=temp_u(ii,jj)
					v_global(i_global,j_global)=temp_v(ii,jj)
					p_global(i_global,j_global)=temp_p(ii,jj)
					xx_global(i_global,j_global)=temp_xx(ii,jj)
					yy_global(i_global,j_global)=temp_yy(ii,jj)
				enddo
				enddo	
				deallocate(temp_rho,temp_u,temp_v,temp_p,temp_xx,temp_yy)		
			endif
100		continue

		open(70,file=trim(filename))
		!open(70,file="initfield.plt")
		write(70,*)'variables=x,y,rho,u,v,p'
		write(70,*)'zone i=',nx_global,'j=',ny_global

		do jj=1,ny_global
		do ii=1,nx_global
		
			write(70,*)xx_global(ii,jj),yy_global(ii,jj),rho_global(ii,jj),u_global(ii,jj),v_global(ii,jj),p_global(ii,jj)
		enddo
		enddo

		close(70)

	else

		allocate(temp_u(i_nn(npx),j_nn(npy)),temp_v(i_nn(npx),j_nn(npy)),temp_p(i_nn(npx),j_nn(npy)), &
				temp_xx(i_nn(npx),j_nn(npy)),temp_yy(i_nn(npx),j_nn(npy)),temp_rho(i_nn(npx),j_nn(npy)))

		do jj=1,j_nn(npy)
		do ii=1,i_nn(npx)
			temp_rho(ii,jj)=rho(ii,jj)		
			temp_u(ii,jj)=u(ii,jj)
			temp_v(ii,jj)=v(ii,jj)
			temp_p(ii,jj)=p(ii,jj)
			temp_xx(ii,jj)=xx(ii,jj)
			temp_yy(ii,jj)=yy(ii,jj)
		enddo
		enddo

		call MPI_SEND(temp_rho,i_nn(npx)*j_nn(npy),OCFD_DATA_TYPE,  &
			0,6,MPI_COMM_WORLD,ierr) 
		call MPI_SEND(temp_u,i_nn(npx)*j_nn(npy),OCFD_DATA_TYPE,  &
			0,1,MPI_COMM_WORLD,ierr) 
		call MPI_SEND(temp_v,i_nn(npx)*j_nn(npy),OCFD_DATA_TYPE,  &
			0,2,MPI_COMM_WORLD,ierr)
		call MPI_SEND(temp_p,i_nn(npx)*j_nn(npy),OCFD_DATA_TYPE,  &
			0,3,MPI_COMM_WORLD,ierr)
		call MPI_SEND(temp_xx,i_nn(npx)*j_nn(npy),OCFD_DATA_TYPE,  &
			0,4,MPI_COMM_WORLD,ierr) 
		call MPI_SEND(temp_yy,i_nn(npx)*j_nn(npy),OCFD_DATA_TYPE,  &
			0,5,MPI_COMM_WORLD,ierr)		

		deallocate(temp_rho,temp_u,temp_v,temp_p,temp_xx,temp_yy)

	endif

	end subroutine


	

