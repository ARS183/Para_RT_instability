module default
	implicit none
	integer,parameter::N1=771,N2=619,NN1=291,NN2=351,NN3=501
	real(kind=8),parameter::gama1=3.8,gama2=3.2,x1=-27,x2=48,y1=-18,y2=32
	real(kind=8),parameter::M0=1.05,Mv1=0.25,Mv2=0.25,gama=1.4,p0=1/gama/M0**2,rho0=1.0,u0=1.0
	real(kind=8)::Ms=M0,ps=p0,lamda1=1.0,lamda2=1.0,hx,hy
end module
program ShockwaveInitial
	use default
	implicit none
	integer::i,j
	real(kind=8)::p(N1,N2),u(N1,N2),v(N1,N2),rho(N1,N2),f1(N1),f2(N2),dx(N1-1),dy(N2-1),x10=-8.0,y10=2.0,x20=-4.0,y20=-2.0,temp1,temp2,ucgrid1,ucgrid2,maxdx,maxdy

	do i=-NN1,NN2		
		if (i<0)then
			f1(i+NN1+1)=x1*(exp(-gama1*i/NN2)-1)/(exp(gama1)-1)
		else
			f1(i+NN1+1)=-x1*(exp(gama1*i/NN2)-1)/(exp(gama1)-1)
		end if
	end	do
	do j=1,NN1+NN2
		dx(j)=f1(j+1)-f1(j)
	end do
	maxdx=maxval(dx)	
	write(*,*)maxdx
	ucgrid1=int((x2+x1)/maxdx)+1
	write(*,*)ucgrid1
	print *,f1(1)
	do i=1,ucgrid1
		f1(i+NN1+NN2+1)=-x1+maxdx*i
	end do
	

	do j=-(NN3-1)/2,(NN3-1)/2
		if (j<0)then
			f2(j+(NN3-1)/2+1)=y1*(exp(-gama2*j/(NN3-1)*2)-1)/(exp(gama2)-1)
		else 
			f2(j+(NN3-1)/2+1)=-y1*(exp(gama2*j/(NN3-1)*2)-1)/(exp(gama2)-1)
		end if
	end	do
	do j=1,NN3-1
		dy(j)=f2(j+1)-f2(j)
	end do
	maxdy=maxval(dy(1:NN3-1))	
	write(*,*)maxdy
	ucgrid2=int((y2+y1)/maxdy)+1
	write(*,*)ucgrid2

	do i=-(N2-1)/2,-(NN3-1)/2-1
		f2(i+(N2-1)/2+1)=y1-(ucgrid2-(i+(N2-1)/2))*maxdy;
	end do   
	do i=-(NN3-1)/2,-1
		f2(i+(N2-1)/2+1)=y1*(exp(-gama2*i/((NN3-1)/2))-1)/(exp(gama2)-1);
	end do    
	do i=0,(NN3-1)/2
		f2(i+(N2-1)/2+1)=-y1*(exp(gama2*i/((NN3-1)/2))-1)/(exp(gama2)-1);
	end do
	do i=(NN3-1)/2+1,(N2-1)/2
		f2(i+(N2-1)/2+1)=-y1+maxdy*(i-(NN3-1)/2);
	end do

	do j=1,N2-1
		dy(j)=f2(j+1)-f2(j)
	end do

	
	Ms=sqrt((2+(gama-1)*M0**2)/(2*gama*M0**2-(gama-1)))
	lamda1=(gama+1)*M0**2/(2+(gama-1)*M0**2)
	lamda2=(2*gama*M0**2-(gama-1))/(gama+1)
	ps=p0*lamda2

	do i=1,N1
		do j=1,N2
			temp1=exp((1-(f1(i)-x10)**2-(f2(j)-y10)**2)/2)
			temp2=exp((1-(f1(i)-x20)**2-(f2(j)-y20)**2)/2)

			if (f1(i)<0)then
				p(i,j)=1/gama/M0**2*(1-(gama-1)/2*Mv1**2*temp1**2)**(gama/(gama-1))+1/gama/M0**2*(1-(gama-1)/2*Mv2**2*temp2**2)**(gama/(gama-1))-p0
				rho(i,j)=(1-(gama-1)/2*Mv1**2*temp1**2)**(1/(gama-1))+(1-(gama-1)/2*Mv2**2*temp2**2)**(1/(gama-1))-rho0
       			u(i,j)=-(f2(j)-y10)*Mv1/M0*temp1+u0-(f2(j)-y20)*Mv2/M0*temp2
				v(i,j)=(f1(i)-x10)*Mv1/M0*temp1+(f1(i)-x20)*Mv2/M0*temp2
			else
	!			p(i,j)=1/gama/M0**2*(1-(gama-1)/2*Mv**2*temp1**2)**(gama/(gama-1))*lamda2
	!			rho(i,j)=(1-(gama-1)/2*Mv**2*temp1**2)**(1/(gama-1))*lamda1
	!     		u(i,j)=(-(f2(j)-y10)*Mv/M0*temp1+u0)/lamda1
	!			v(i,j)=(f1(i)-x10)*Mv/M0*temp1

				p(i,j)=(1/gama/M0**2*(1-(gama-1)/2*Mv1**2*temp1**2)**(gama/(gama-1))+1/gama/M0**2*(1-(gama-1)/2*Mv2**2*temp2**2)**(gama/(gama-1))-p0)*lamda2
				rho(i,j)=((1-(gama-1)/2*Mv1**2*temp1**2)**(1/(gama-1))+(1-(gama-1)/2*Mv2**2*temp2**2)**(1/(gama-1))-rho0)*lamda1
	      		u(i,j)=(-(f2(j)-y10)*Mv1/M0*temp1+u0-(f2(j)-y20)*Mv2/M0*temp2)/lamda1
				v(i,j)=(f1(i)-x10)*Mv1/M0*temp1+(f1(i)-x20)*Mv2/M0*temp2

			end if			
	   	end do
	end do
		
	p(1,:)=p0
	rho(1,:)=rho0
	u(1,:)=u0
	v(1,:)=0

	open(unit=5,file="u.dat")
	write(5,*)u
	open(unit=10,file="v.dat")
	write(10,*)v
	open(unit=15,file="p.dat")
	write(15,*)p
	open(unit=20,file="rho.dat")
	write(20,*)rho
	hx=1.d0
	hy=1.d0
	call output_driven(rho(:,:),u(:,:),v(:,:),p(:,:),f1(:),f2(:),N1,N2,0)
end


subroutine output_driven(rho,u,v,p,f1,f2,ni,nj,itp)
      implicit double precision(a-h,o-z)
	dimension rho(ni,nj),u(ni,nj),v(ni,nj),p(ni,nj),f1(ni),f2(nj)
	character*6 name2
	character*20 fname2
    
      write(name2,'(I6.6)') int(itp)
      fname2='output_'//'_'//name2//'.plt' 
      open(2,file=fname2,status='unknown')
      write(2,*) 'variables=x,y,rho,u,v,p'
	write(2,*) 'zone i=',ni,'j=',nj
      
	do j=1,nj	
	do  i=1,ni
	x=f1(i)
	y=f2(j)
	write(2,*)x,y,rho(i,j),u(i,j),v(i,j),p(i,j)
	end do
	end do
    close(2)
	end subroutine