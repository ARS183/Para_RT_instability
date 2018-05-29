module default
	implicit none
	integer,parameter::N1=771,N2=619,NN1=291,NN2=351,NN3=501,M=200000
	real(kind=8),parameter::gama1=3.8,gama2=3.2,x1=-27,x2=48,y1=-18,y2=32
	real(kind=8),parameter::h1=0.1,h2=0.1,dt=5E-4,sigma=0.54395,Re=500,M0=1.05,Cu=110.4,  &
	gama=1.4,p0=1/gama/M0**2,rho0=1.0,Pr=0.72,u0=1.0
	real(kind=8)::ex(N1,N2)=1.0,ny(N1,N2)=1.0,Jac(N1,N2)=1.0,  &
	Jac1(N1,N2)=1.0,lamda1=1.0,lamda2=1.0
	type FLUX
		real(kind=8)::Fir(N1,N2)
		real(kind=8)::Sec(N1,N2)
		real(kind=8)::Thir(N1,N2)
		real(kind=8)::Four(N1,N2)
	end type FLUX
end module
program ShockWave
	use default
	implicit none
	interface
		function TimePace(u,v,p,rho)
			use default
			real(kind=8)::p(N1,N2),u(N1,N2),v(N1,N2),rho(N1,N2)
			type(FLUX)::TimePace
		end function
	end interface

	integer::i,j
	real(kind=8)::p(N1,N2),u(N1,N2),v(N1,N2),rho(N1,N2)
	type(FLUX)::Tmp	
	character(len=8)::filename1,filename2,filename3,filename4

	filename1="u000.dat"
	filename2="v000.dat"
	filename3="p000.dat"
	filename4="r000.dat"

	do j=1,N2
		do i=-NN1,N1-NN1-1		
			if (i<0)then
				ex(i+NN1+1,j)=(exp(gama1)-1)/(-x1/NN2*gama1)*exp(i*gama1/NN2)*h1
			else if(i<=NN2)then
				ex(i+NN1+1,j)=(exp(gama1)-1)/(-x1/NN2*gama1)*exp(-i*gama1/NN2)*h1
			else
				ex(i+NN1+1,j)=ex(NN1+NN2+1,1)
			end if
		end	do
	end do

	do i=1,N1
		do j=-(N2-1)/2,-(NN3-1)/2-1
			ny(i,j+(N2-1)/2+1)=(exp(gama2)-1)/(-y1/(NN3-1)*2*gama2)*  &
			exp(-(NN3-1)/2*gama2/(NN3-1)*2)*h2
		end do
		do j=-(NN3-1)/2,-1
			ny(i,j+(N2-1)/2+1)=(exp(gama2)-1)/(-y1/(NN3-1)*2*gama2)*exp(j*gama2/(NN3-1)*2)*h2
		end do
		do j=0,(NN3-1)/2
			ny(i,j+(N2-1)/2+1)=(exp(gama2)-1)/(-y1/(NN3-1)*2*gama2)*exp(-j*gama2/(NN3-1)*2)*h2
		end do
		do j=(NN3-1)/2+1,(N2-1)/2
			ny(i,j+(N2-1)/2+1)=ny(1,(NN3-1)/2+(N2-1)/2+1)			
		end	do
	end do

!	open(unit=15,file="ex.dat")
!	write(15,*)ex
!	open(unit=20,file="ny.dat")
!	write(20,*)ny
!	pause

	Jac=ex*ny
	Jac1=ex/ny

	lamda1=(gama+1)*M0**2/(2+(gama-1)*M0**2)
	lamda2=(2*gama*M0**2-(gama-1))/(gama+1)

	open(unit=15,file="u.dat",status="old")
	read(15,*)u
	open(unit=20,file="v.dat",status="old")
	read(20,*)v
	open(unit=25,file="p.dat",status="old")
	read(25,*)p
	open(unit=30,file="rho.dat",status="old")
	read(30,*)rho
	close(15)
	close(20)
	close(25)
	close(30)



	do i=1,M
		write(*,*)i
		Tmp=TimePace(u,v,p,rho)
		rho=Tmp%Fir
		u=Tmp%Sec
		v=Tmp%Thir
		p=Tmp%four
		p(1,:)=p0
		rho(1,:)=rho0
		u(1,:)=u0
		v(1,:)=0
		if (mod(i,1000)==0)then

			j=i/1000
			filename1(2:2)=achar(j/100+48)
			filename2(2:2)=achar(j/100+48)
			filename3(2:2)=achar(j/100+48)
			filename4(2:2)=achar(j/100+48)

			filename1(3:3)=achar(mod(j/10,10)+48)
			filename2(3:3)=achar(mod(j/10,10)+48)
			filename3(3:3)=achar(mod(j/10,10)+48)
			filename4(3:3)=achar(mod(j/10,10)+48)

			filename1(4:4)=achar(j-j/100*100-mod(j/10,10)*10+48)
			filename2(4:4)=achar(j-j/100*100-mod(j/10,10)*10+48)
			filename3(4:4)=achar(j-j/100*100-mod(j/10,10)*10+48)
			filename4(4:4)=achar(j-j/100*100-mod(j/10,10)*10+48)

			open(unit=1,file=filename1)
			open(unit=2,file=filename2)
			open(unit=3,file=filename3)
			open(unit=4,file=filename4)

			write(1,*)u
			write(2,*)v
			write(3,*)p

			write(4,*)rho
			close(1)
			close(2)
			close(3)
			close(4)
		end if
	end do

!	open(unit=35,file="u.dat",status="replace")
!	write(35,*)u
!	open(unit=40,file="v.dat",status="replace")
!	write(40,*)v
!	open(unit=45,file="P.dat",status="replace")
!	write(45,*)p
!	open(unit=50,file="rho.dat",status="replace")
!	write(50,*)rho
	
	
end


subroutine output_driven(rho,u,v,f1,f2,ni,nj,hx,hy,itp)
      implicit double precision(a-h,o-z)
	dimension rho(ni,nj),u(ni,nj),v(ni,nj),f1(ni,nj),f2(ni,nj)
	character*6 name2
	character*20 fname2
	integer:: itp
    
      write(name2,'(I6.6)') int(itp)
      fname2='output_'//'_'//name2//'.plt' 
      open(2,file=fname2,status='unknown')
      write(2,*) 'variables=x,y,rho,u,v'
	write(2,*) 'zone i=',ni,'j=',nj
      
	do j=1,nj	
	do  i=1,ni
	x=f1(i,j)
	y=f2(i,j)
	write(2,*)x,y,rho(i,j),u(i,j),v(i,j)
	end do
	end do
    close(2)
	end subroutine
function Chas(a,b,c,d,N)
	implicit none
	integer::N,j
	real(kind=8)::a(N),b(N),c(N),d(N)
	real(kind=8)::Chas(N)
    do j=2,N
        a(j)=a(j)/b(j-1)
        b(j)=b(j)-c(j-1)*a(j)
        d(j)=d(j)-a(j)*d(j-1)
    end	do
    d(N)=d(N)/b(N)
    do j=N-1,1,-1
        d(j)=(d(j)-c(j)*d(j+1))/b(j)
    end	do
	Chas=d
	return
end
function CirChas(a,b,c,d,N)
	implicit none
	integer::N
	integer::j
	real(kind=8)::a(N),b(N),c(N),d(N)
	real(kind=8)::p(N),q(N),r(N-1),t(N-1),uu(N),vv(N),s
	real(kind=8)::CirChas(N)

	p(1)=b(1)
	q(1)=c(1)/p(1)
	r(1)=a(1)/b(1)
    do j=2,N-2
        p(j)=b(j)-a(j)*q(j-1)
        q(j)=c(j)/p(j)
		r(j)=-a(j)*r(j-1)/p(j)
    end do
    p(N-1)=b(N-1)-a(N-1)*q(N-2)
    r(N-1)=(c(N-1)-a(N-1)*r(N-2))/p(N-1)
	p(N)=b(N)-a(N)*r(N-1)
	s=c(N)/p(N)

	t(N-1)=r(N-1)
    do j=N-2,1,-1 
        t(j)=r(j)-q(j)*t(j+1)
    end do
	
	uu(1)=d(1)/p(1)
	do j=2,N
        uu(j)=(d(j)-a(j)*uu(j-1))/p(j)
    end do

	vv(N)=uu(N)
	vv(N-1)=uu(N-1)
	do j=N-2,1,-1
        vv(j)=uu(j)-q(j)*vv(j+1)
    end do

	CirChas(N)=(vv(N)-s*vv(1))/(1-t(1)*s)
	do j=N-1,1,-1
        CirChas(j)=vv(j)-t(j)*CirChas(N)
	end do

	return
end
function CCUpos1(h,u,N)
	use default
	implicit none
	interface
		function CirChas(a,b,c,d,N)
			integer::N
			real(kind=8)::a(N),b(N),c(N),d(N)
			real(kind=8)::CirChas(N)
		end function
	end interface

	integer::i,j,N
	real(kind=8)::h
	real(kind=8)::u(N)
	real(kind=8)::a(N),b(N),c(N),d(N),S(N)
	real(kind=8)::CCUpos1(N)

	a=1.0/12.0
	b=5.0/6.0
	c=1.0/12.0

    d(1)=(u(N)-2.0*u(1)+u(2))/h**2
    do j=2,N-1
        d(j)=(u(j-1)-2.0*u(j)+u(j+1))/h**2
    end do
    d(N)=(u(N-1)-2.0*u(N)+u(1))/h**2
    S=CirChas(a,b,c,d,N)
    
	a=1.0+sigma
	b=2.0+sigma
	c=0.0
   
    d(1)=(-3.0*sigma*u(N-1)+(-42.0-2.0*sigma)*u(N)+(48.0-12.0*sigma)*u(1)+(-6.0+18.0*sigma)*u(2)-sigma*u(3))/12.0/h+h*(1.0-sigma)*S(1)
    d(2)=(-3.0*sigma*u(N)+(-42.0-2.0*sigma)*u(1)+(48.0-12.0*sigma)*u(2)+(-6.0+18.0*sigma)*u(3)-sigma*u(4))/12.0/h+h*(1.0-sigma)*S(2)
    do j=3,N-2
        d(j)=(-3.0*sigma*u(j-2)+(-42.0-2.0*sigma)*u(j-1)+(48.0-12.0*sigma)*u(j)+(-6.0+18.0*sigma)*u(j+1)-sigma*u(j+2))/12.0/h+h*(1.0-sigma)*S(j)
    end	do
    d(N-1)=(-3.0*sigma*u(N-3)+(-42.0-2.0*sigma)*u(N-2)+(48.0-12.0*sigma)*u(N-1)+(-6.0+18.0*sigma)*u(N)-sigma*u(1))/12.0/h+h*(1.0-sigma)*S(N-1)
    d(N)=(-3.0*sigma*u(N-2)+(-42.0-2.0*sigma)*u(N-1)+(48.0-12.0*sigma)*u(N)+(-6.0+18.0*sigma)*u(1)-sigma*u(2))/12.0/h+h*(1.0-sigma)*S(N)

    CCUpos1=CirChas(a,b,c,d,N)
	return
end
function CCUpos2(h,u,N)   
	use default
	implicit none
	interface
		function Chas(a,b,c,d,N)
			integer::N
			real(kind=8)::a(N),b(N),c(N),d(N)
			real(kind=8)::Chas(N)
		end function
	end interface

	integer::i,j,N
	real(kind=8)::h
	real(kind=8)::u(N)
	real(kind=8)::a(N-2),b(N-2),c(N-2),d(N-2),S(N)
	real(kind=8)::CCUpos2(N)

    a=1.0/12.0
    b=5.0/6.0
    c=1.0/12.0
	b(1)=1.0/12.0+5.0/6.0
	b(N-2)=1.0/12.0+5.0/6.0
	
    do j=2,N-1
        d(j-1)=(u(j-1)-2*u(j)+u(j+1))/h**2
    end do
    S(2:N-1)=Chas(a,b,c,d,N-2)
	
    CCUpos2(1)=(-3*u(1)+4*u(2)-u(3))/2/h
    CCUpos2(2)=((-5*u(1)+4*u(2)+u(3))/6.0/h-CCUpos2(1)/3)*1.5
    do j=3,N-2
        CCUpos2(j)=((-3*sigma*u(j-2)+(-42-2*sigma)*u(j-1)+(48-12*sigma)*u(j)+(-6+18*sigma)*u(j+1)-sigma*u(j+2))/12/h+h*(1-sigma)*S(j)-(1+sigma)*CCUpos2(j-1))/(2+sigma)
    end do
    CCUpos2(N-1)=((-5*u(N-2)+4*u(N-1)+u(N))/6.0/h-CCUpos2(N-2)/3)*1.5
    CCUpos2(N)=(3*u(N)-4*u(N-1)+u(N-2))/2/h
	
end
function CCUneg1(h,u,N)
	use default
	implicit none
	interface
		function CirChas(a,b,c,d,N)
			integer::N
			real(kind=8)::a(N),b(N),c(N),d(N)
			real(kind=8)::CirChas(N)
		end function
	end interface

	integer::i,j,N
	real(kind=8)::h
	real(kind=8)::u(N)
	real(kind=8)::a(N),b(N),c(N),d(N),S(N),FF(N)
	real(kind=8)::CCUneg1(N)

	a=1.0/12.0
	b=5.0/6.0
	c=1.0/12.0

    d(1)=(u(N)-2.0*u(1)+u(2))/h**2
    do j=2,N-1
        d(j)=(u(j-1)-2.0*u(j)+u(j+1))/h**2
    end do
    d(N)=(u(N-1)-2.0*u(N)+u(1))/h**2
    S=CirChas(a,b,c,d,N-2)
	
	do i=1,N
		a(i)=1.0+sigma
		b(i)=2.0+sigma
		c(i)=0.0
	end do

	d(1)=(sigma*u(N-2)+(6.0-18.0*sigma)*u(N-1)+(-48.0+12.0*sigma)*u(N)+(42.0+2.0*sigma)*u(1)+3.0*sigma*u(2))/12.0/h-h*(1.0-sigma)*S(N)
    d(2)=(sigma*u(N-3)+(6.0-18.0*sigma)*u(N-2)+(-48.0+12.0*sigma)*u(N-1)+(42.0+2.0*sigma)*u(N)+3.0*sigma*u(1))/12.0/h-h*(1.0-sigma)*S(N-1)
    do j=3,N-2
        d(j)=(sigma*u(N-1-j)+(6.0-18.0*sigma)*u(N-j)+(-48.0+12.0*sigma)*u(N+1-j)+(42.0+2.0*sigma)*u(N+2-j)+3.0*sigma*u(N+3-j))/12.0/h-h*(1.0-sigma)*S(N+1-j)
    end do
    d(N-1)=(sigma*u(N)+(6.0-18.0*sigma)*u(1)+(-48.0+12.0*sigma)*u(2)+(42.0+2.0*sigma)*u(3)+3.0*sigma*u(4))/12.0/h-h*(1.0-sigma)*S(2)
    d(N)=(sigma*u(N-1)+(6.0-18.0*sigma)*u(N)+(-48.0+12.0*sigma)*u(1)+(42.0+2.0*sigma)*u(2)+3.0*sigma*u(3))/12.0/h-h*(1.0-sigma)*S(1)
        
    FF=CirChas(a,b,c,d,N)    
    do j=1,N
        CCUneg1(j)=FF(N+1-j)
    end	do
	return
end
function CCUneg2(h,u,N)
	use default
	implicit none
	interface
		function Chas(a,b,c,d,N)
			integer::N
			real(kind=8)::a(N),b(N),c(N),d(N)
			real(kind=8)::Chas(N)
		end function
	end interface

	integer::i,j,N
	real(kind=8)::h
	real(kind=8)::u(N)
	real(kind=8)::a(N-2),b(N-2),c(N-2),d(N-2),S(N)
	real(kind=8)::CCUneg2(N)

    a=1.0/12.0
    b=5.0/6.0
    c=1.0/12.0
	b(1)=1.0/12.0+5.0/6.0
	b(N-2)=1.0/12.0+5.0/6.0
	
    do j=2,N-1
        d(j-1)=(u(j-1)-2*u(j)+u(j+1))/h**2
    end do
    S(2:N-1)=Chas(a,b,c,d,N-2)
	
    CCUneg2(N)=(3*u(N)-4*u(N-1)+u(N-2))/2/h
    CCUneg2(N-1)=((5*u(N)-4*u(N-1)-u(N-2))/6.0/h-CCUneg2(N)/3)*1.5
    do j=N-2,3,-1
        CCUneg2(j)=((sigma*u(j-2)+(6-18*sigma)*u(j-1)+(-48+12*sigma)*u(j)+(42+2*sigma)*u(j+1)+3*sigma*u(j+2))/12/h-h*(1-sigma)*S(j)-(1+sigma)*CCUneg2(j+1))/(2+sigma)
    end do
    CCUneg2(2)=((5*u(3)-4*u(2)-u(1))/6.0/h-CCUneg2(3)/3)*1.5
    CCUneg2(1)=(-3*u(1)+4*u(2)-u(3))/2/h
	
    return
end

function VanLeerE(u,v,p,rho,cs,Ep)
	use default
	implicit none
	interface
		function CCUpos2(h,u,N)
			integer::N
			real(kind=8)::h
			real(kind=8)::u(N)
			real(kind=8)::CCUpos2(N)
		end function CCUpos2
		function CCUneg2(h,u,N)
			integer::N
			real(kind=8)::h
			real(kind=8)::u(N)
			real(kind=8)::CCUneg2(N)
		end function CCUneg2
	end interface
	real(kind=8)::p(N1,N2),u(N1,N2),v(N1,N2),rho(N1,N2),cs(N1,N2),Me(N1,N2)
	real(kind=8)::unit1
	integer::i,j
	type(FLUX)::Ep,Epos,Eneg,VanLeerE

	Me=u/cs

	do i=1,N1
		do j=1,N2
			if(Me(i,j)>=1)then
				Epos%Fir(i,j)=Ep%Fir(i,j)
				Epos%Sec(i,j)=Ep%Sec(i,j)
				Epos%Thir(i,j)=Ep%Thir(i,j)
				Epos%Four(i,j)=Ep%Four(i,j)

				Eneg%Fir(i,j)=0
				Eneg%Sec(i,j)=0
				Eneg%Thir(i,j)=0
				Eneg%Four(i,j)=0
			else if(Me(i,j)<=-1)then
				Epos%Fir(i,j)=0
				Epos%Sec(i,j)=0
				Epos%Thir(i,j)=0
				Epos%Four(i,j)=0

				Eneg%Fir(i,j)=Ep%Fir(i,j)
				Eneg%Sec(i,j)=Ep%Sec(i,j)
				Eneg%Thir(i,j)=Ep%Thir(i,j)
				Eneg%Four(i,j)=Ep%Four(i,j)
			else
				unit1=(-(gama-1)*u(i,j)**2+2*(gama-1)*u(i,j)*cs(i,j)+2.0*cs(i,j)**2)/(gama**2-1)+(u(i,j)**2+v(i,j)**2)/2.0

				Epos%Fir(i,j)=rho(i,j)*cs(i,j)*(Me(i,j)+1)**2/ny(i,j)/4.0
				Epos%Sec(i,j)=Epos%Fir(i,j)*((-u(i,j)+2*cs(i,j))/gama+u(i,j))
				Epos%Thir(i,j)=Epos%Fir(i,j)*v(i,j)
				Epos%Four(i,j)=Epos%Fir(i,j)*unit1				
			
				Eneg%Fir(i,j)=-rho(i,j)*cs(i,j)*(Me(i,j)-1)**2/ny(i,j)/4.0
				Eneg%Sec(i,j)=Eneg%Fir(i,j)*((-u(i,j)-2*cs(i,j))/gama+u(i,j))
				Eneg%Thir(i,j)=Eneg%Fir(i,j)*v(i,j)
				Eneg%Four(i,j)=Eneg%Fir(i,j)*(unit1-4*u(i,j)*cs(i,j)/(gama+1))
			end if
		end do
	end do
	do j=1,N2	
		VanLeerE%Fir(:,j)=CCUpos2(h1,Epos%Fir(:,j),N1)+CCUneg2(h1,Eneg%Fir(:,j),N1)
		VanLeerE%Sec(:,j)=CCUpos2(h1,Epos%Sec(:,j),N1)+CCUneg2(h1,Eneg%Sec(:,j),N1)
		VanLeerE%Thir(:,j)=CCUpos2(h1,Epos%Thir(:,j),N1)+CCUneg2(h1,Eneg%Thir(:,j),N1)
		VanLeerE%Four(:,j)=CCUpos2(h1,Epos%Four(:,j),N1)+CCUneg2(h1,Eneg%Four(:,j),N1)
	end do
	return
end
function VanLeerF(u,v,p,rho,cs,Fp)
	use default
	implicit none
		interface
		function CCUpos2(h,u,N)
			integer::N
			real(kind=8)::h
			real(kind=8)::u(N)
			real(kind=8)::CCUpos2(N)
		end function CCUpos2
		function CCUneg2(h,u,N)
			integer::N
			real(kind=8)::h
			real(kind=8)::u(N)
			real(kind=8)::CCUneg2(N)
		end function CCUneg2
	end interface
	real(kind=8)::p(N1,N2),u(N1,N2),v(N1,N2),rho(N1,N2),cs(N1,N2),Mn(N1,N2)
	real(kind=8)::unit1
	integer::i,j
	type(FLUX)::Fp,Fpos,Fneg,VanLeerF

	Mn=v/cs

	do i=1,N1
		do j=1,N2
			if(Mn(i,j)>=1)then
				Fpos%Fir(i,j)=Fp%Fir(i,j)
				Fpos%Sec(i,j)=Fp%Sec(i,j)
				Fpos%Thir(i,j)=Fp%Thir(i,j)
				Fpos%Four(i,j)=Fp%Four(i,j)

				Fneg%Fir(i,j)=0
				Fneg%Sec(i,j)=0
				Fneg%Thir(i,j)=0
				Fneg%Four(i,j)=0
			else if(Mn(i,j)<=-1)then
				Fpos%Fir(i,j)=0
				Fpos%Sec(i,j)=0
				Fpos%Thir(i,j)=0
				Fpos%Four(i,j)=0

				Fneg%Fir(i,j)=Fp%Fir(i,j)
				Fneg%Sec(i,j)=Fp%Sec(i,j)
				Fneg%Thir(i,j)=Fp%Thir(i,j)
				Fneg%Four(i,j)=Fp%Four(i,j)
			else
				unit1=(-(gama-1)*v(i,j)**2+2*(gama-1)*v(i,j)*cs(i,j)+2.0*cs(i,j)**2)/(gama**2-1)+(u(i,j)**2+v(i,j)**2)/2.0

				Fpos%Fir(i,j)=rho(i,j)*cs(i,j)*(Mn(i,j)+1)**2/ex(i,j)/4.0
				Fpos%Sec(i,j)=Fpos%Fir(i,j)*u(i,j)
				Fpos%Thir(i,j)=Fpos%Fir(i,j)*((-v(i,j)+2*cs(i,j))/gama+v(i,j))
				Fpos%Four(i,j)=Fpos%Fir(i,j)*unit1				
			
				Fneg%Fir(i,j)=-rho(i,j)*cs(i,j)*(Mn(i,j)-1)**2/ex(i,j)/4.0
				Fneg%Sec(i,j)=Fneg%Fir(i,j)*u(i,j)
				Fneg%Thir(i,j)=Fneg%Fir(i,j)*((-v(i,j)-2*cs(i,j))/gama+v(i,j))
				Fneg%Four(i,j)=Fneg%Fir(i,j)*(unit1-4*v(i,j)*cs(i,j)/(gama+1))
			end if
		end do
	end do

	do i=1,N1
		VanLeerF%Fir(i,:)=CCUpos2(h2,Fpos%Fir(i,:),N2)+CCUneg2(h2,Fneg%Fir(i,:),N2)
		VanLeerF%Sec(i,:)=CCUpos2(h2,Fpos%Sec(i,:),N2)+CCUneg2(h2,Fneg%Sec(i,:),N2)
		VanLeerF%Thir(i,:)=CCUpos2(h2,Fpos%Thir(i,:),N2)+CCUneg2(h2,Fneg%Thir(i,:),N2)
		VanLeerF%Four(i,:)=CCUpos2(h2,Fpos%Four(i,:),N2)+CCUneg2(h2,Fneg%Four(i,:),N2)
	end do
	return
end
function SCD61(h,N,u) 
	implicit none
	interface
		function CirChas(a,b,c,d,N)
			integer::N
			real(kind=8)::a(N),b(N),c(N),d(N)
			real(kind=8)::CirChas(N)
		end function
	end interface

	integer::N
	real(kind=8)::h
	real(kind=8)::u(N)
	real(kind=8)::SCD61(N)	
	real(kind=8)::a(N),b(N),c(N),d(N),S(N)
	integer::i,j

	a=1.0/3.0
	b=1.0
	c=1.0/3.0


	d(1)=(u(3)+28*u(2)-28*u(N)-u(N-1))/36.0/h
	d(2)=(u(4)+28*u(3)-28*u(1)-u(N))/36.0/h
    do j=3,N-2
        d(j)=(u(j+2)+28*u(j+1)-28*u(j-1)-u(j-2))/36.0/h
    end do
	d(N-1)=(u(1)+28*u(N)-28*u(N-2)-u(N-3))/36.0/h
    d(N)=(u(2)+28*u(1)-28*u(N-1)-u(N-2))/36.0/h
    SCD61=CirChas(a,b,c,d,N)
end
function SCD62(h,N,u) 
	implicit none
	interface
		function Chas(a,b,c,d,N)
			integer::N
			real(kind=8)::a(N),b(N),c(N),d(N)
			real(kind=8)::Chas(N)
		end function
	end interface

	integer::N
	real(kind=8)::h
	real(kind=8)::u(N)
	real(kind=8)::SCD62(N)	
	real(kind=8)::a(N),b(N),c(N),d(N),S(N)
	integer::i,j

	a=1.0/3.0
	b=1.0
	c=1.0/3.0

	b(1)=1.0/3.0
	b(2)=2.0/3.0
	b(N-1)=2.0/3.0
	b(N)=1.0/3.0

	a(2)=1.0/4.0
	a(N-1)=1.0/4.0
	a(N)=2.0/3.0

	c(1)=2.0/3.0
	c(2)=1.0/4.0
	c(N-1)=1.0/4.0

	d(1)=(-5*u(1)+4*u(2)+u(3))/6.0/h
	d(2)=(u(3)-u(1))*0.75/h
    do j=3,N-2
        d(j)=(u(j+2)+28*u(j+1)-28*u(j-1)-u(j-2))/36.0/h
    end do
	d(N-1)=(u(N)-u(N-2))*0.75/h
    d(N)=(5*u(N)-4*u(N-1)-u(N-2))/6.0/h
    SCD62=Chas(a,b,c,d,N)
end

function Ru(u,v,p,rho)
	use default
	implicit none
	interface
		function SCD61(h,N,u)
			real(kind=8)::h
			real(kind=8)::u(N)
			real(kind=8)::SCD61(N)
		end function SCD61
		function SCD62(h,N,u)
			real(kind=8)::h
			real(kind=8)::u(N)
			real(kind=8)::SCD62(N)
		end function SCD62
		function VanleerE(u,v,p,rho,cs,Ep)
			use default
			real(kind=8)::p(N1,N2),u(N1,N2),v(N1,N2),rho(N1,N2),cs(N1,N2)
			type(FLUX)::Ep,VanleerE
		end	function VanLeerE
		function VanleerF(u,v,p,rho,cs,Fp)
			use default
			real(kind=8)::p(N1,N2),u(N1,N2),v(N1,N2),rho(N1,N2),cs(N1,N2)
			type(FLUX)::Fp,VanleerF	
		end	function VanLeerF				
	end interface
	
	integer::i,j
	real(kind=8)::p(N1,N2),u(N1,N2),v(N1,N2),rho(N1,N2),T(N1,N2),eng(N1,N2),miu(N1,N2),k(N1,N2),cs(N1,N2)
	real(kind=8)::uu(N1,N2),vv(N1,N2),ue(N1,N2),un(N1,N2),ve(N1,N2),vn(N1,N2),Te(N1,N2),Tn(N1,N2)

	type(FLUX)::Ep,Fp,Ee,Fn,Eve,Evn,Fve,Fvn,Evee,Evne,Fven,Fvnn
	type(FLUX)::Ru

	eng=p/(gama-1.0)+rho/2.0*(u**2+v**2)
	T=p/rho/(gama-1)
    miu=T**1.5*(1.0+Cu/288.15)/(T+Cu/288.15)
    k=miu*gama/Pr
	cs=sqrt(gama*p/rho)

	Ep%Fir=rho*u/ny
	Ep%Sec=(rho*u**2+p)/ny
	Ep%Thir=(rho*u*v)/ny
	Ep%Four=(eng+p)*u/ny

	Fp%Fir=rho*v/ex
	Fp%Sec=(rho*u*v)/ex
	Fp%Thir=(rho*v**2+p)/ex
	Fp%Four=(eng+p)*v/ex

	Ee=VanLeerE(u,v,p,rho,cs,Ep)
	Fn=VanleerF(u,v,p,rho,cs,Fp)

	do j=1,N2
		ue(:,j)=SCD62(h1,N1,u(:,j))
	end do
	do i=1,N1
		un(i,:)=SCD62(h2,N2,u(i,:))
	end do
	do j=1,N2
		ve(:,j)=SCD62(h1,N1,v(:,j))
	end do
	do i=1,N1
		vn(i,:)=SCD62(h2,N2,v(i,:))
	end do
	do j=1,N2
		Te(:,j)=SCD62(h1,N1,T(:,j))
	end do
	do i=1,N1
		Tn(i,:)=SCD62(h2,N2,T(i,:))
	end do

	Eve%Fir=0
	Eve%Sec=Jac1*4/3*miu*ue 
	Eve%Thir=Jac1*miu*ve
	Eve%Four=Jac1*(4/3*miu*u*ue+miu*v*ve+k*Te)

	Evn%Fir=0 
	Evn%Sec=-2/3*miu*vn 
	Evn%Thir=miu*un 
	Evn%Four=-2/3*miu*u*vn+miu*v*un

	Fve%Fir=0 
	Fve%Sec=miu*ve
	Fve%Thir=-2/3*miu*ue
	Fve%Four=miu*u*ve-2/3*miu*v*ue

	Fvn%Fir=0
	Fvn%Sec=miu*un/Jac1
	Fvn%Thir=4/3*miu*vn/Jac1
	Fvn%Four=(miu*u*un+4/3*miu*v*vn+k*Tn)/Jac1

	do j=1,N2
		Evee%Fir(:,j)=SCD62(h1,N1,Eve%Fir(:,j))
		Evee%Sec(:,j)=SCD62(h1,N1,Eve%Sec(:,j))
		Evee%Thir(:,j)=SCD62(h1,N1,Eve%Thir(:,j))
		Evee%Four(:,j)=SCD62(h1,N1,Eve%Four(:,j))
	end do
	do j=1,N2
		Evne%Fir(:,j)=SCD62(h1,N1,Evn%Fir(:,j))
		Evne%Sec(:,j)=SCD62(h1,N1,Evn%Sec(:,j))
		Evne%Thir(:,j)=SCD62(h1,N1,Evn%Fir(:,j))
		Evne%Four(:,j)=SCD62(h1,N1,Evn%Four(:,j))
	end do
	do i=1,N1
		Fven%Fir(i,:)=SCD62(h2,N2,Fve%Fir(i,:))
		Fven%Sec(i,:)=SCD62(h2,N2,Fve%Sec(i,:))
		Fven%Thir(i,:)=SCD62(h2,N2,Fve%Thir(i,:))
		Fven%Four(i,:)=SCD62(h2,N2,Fve%Four(i,:))
	end do
	do i=1,N1
		Fvnn%Fir(i,:)=SCD62(h2,N2,Fvn%Fir(i,:))
		Fvnn%Sec(i,:)=SCD62(h2,N2,Fvn%Sec(i,:))
		Fvnn%Thir(i,:)=SCD62(h2,N2,Fvn%Thir(i,:))
		Fvnn%Four(i,:)=SCD62(h2,N2,Fvn%Four(i,:))
	end do

	RU%Fir=-Ee%Fir-Fn%Fir +(Evee%Fir+Evne%Fir+Fven%Fir+Fvnn%Fir)/Re
	RU%Sec=-Ee%Sec-Fn%Sec +(Evee%Sec+Evne%Sec+Fven%Sec+Fvnn%Sec)/Re
	RU%Thir=-Ee%Thir-Fn%Thir +(Evee%Thir+Evne%Thir+Fven%Thir+Fvnn%Thir)/Re
	RU%Four=-Ee%Four-Fn%Four +(Evee%Four+Evne%Four+Fven%Four+Fvnn%Four)/Re

	return
end

function TimePace(u,v,p,rho)
	use default
	implicit none
	interface
		function Ru(u,v,p,rho)
			use default
			real(kind=8)::p(N1,N2),u(N1,N2),v(N1,N2),rho(N1,N2)
			type(FLUX)Ru
		end function
	end interface
	integer::i,j
	real(kind=8)::p(N1,N2),u(N1,N2),v(N1,N2),rho(N1,N2),eng(N1,N2)
	type(FLUX)::Up,Up1,Up2,Up3,Rn1,Rn2,Rn3,Rn,Temp
	type(FLUX)::TimePace
	real(kind=8)::a11,a21,a22,a31,a32,a33,b1,b2,b3,c1,c2,c3,error

	a11= 0.377847764031163D0
	a21= 0.385232756462588D0
	a22= 0.461548399939329D0
	a31= 0.675724855841358D0
	a32=-0.061710969841169D0
	a33= 0.241480233100410D0
	b1=  0.750869573741408D0
	b2= -0.362218781852651D0
	b3=  0.611349208111243D0
	c1=  0.257820901066211D0
	c2=  0.434296446908075D0
	c3=  0.758519768667167D0

	eng=p/(gama-1.0)+rho/2.0*(u**2+v**2)

	Up%Fir=rho/Jac
	Up%Sec=rho*u/Jac
	Up%Thir=rho*v/Jac
	Up%Four=eng/Jac

	
!=========================F1==================================
	error=1
	Up1%Fir=Up%Fir
	Up1%Sec=Up%Sec
	Up1%Thir=Up%Thir
	Up1%Four=Up%Four

	do while(abs(error)>10e-7)
	error=0
	Temp%Fir=Up1%Fir
	Temp%Sec=Up1%Sec
	Temp%Thir=Up1%Thir
	Temp%Four=Up1%Four

	Rn=RU(u,v,p,rho)
    Up1%Fir=Up%Fir+dt*a11*Rn%Fir
    Up1%Sec=Up%Sec+dt*a11*Rn%Sec
    Up1%Thir=Up%Thir+dt*a11*Rn%Thir
    Up1%Four=Up%Four+dt*a11*Rn%Four

    rho=Jac*Up1%Fir
    u=Jac*Up1%Sec/rho
    v=Jac*Up1%Thir/rho
    p=(gama-1.0)*(Jac*Up1%Four-rho/2.0*(u**2+v**2))
	p(1,:)=p0
	rho(1,:)=rho0
	u(1,:)=u0
	v(1,:)=0

	Up1%Fir=rho/Jac
	Up1%Sec=rho*u/Jac
	Up1%Thir=rho*v/Jac
	Up1%Four=eng/Jac

	do i=1,N1
	do j=1,N2
	 if (abs(Up1%Fir(i,j)-Temp%Fir(i,j))>error) then
	 error=abs(Up1%Fir(i,j)-Temp%Fir(i,j))
	 end if
	if (abs(Up1%Sec(i,j)-Temp%Sec(i,j))>error) then
	 error=abs(Up1%Sec(i,j)-Temp%Sec(i,j))
	 end if
	if (abs(Up1%Thir(i,j)-Temp%Thir(i,j))>error) then
	 error=abs(Up1%Thir(i,j)-Temp%Thir(i,j))
	 end if
	if (abs(Up1%Four(i,j)-Temp%Four(i,j))>error) then
	 error=abs(Up1%Four(i,j)-Temp%Four(i,j))
	 end if
	end do
	end do
	print *,error
	end do


!=======================F2==================================
	error=1
	Up2%Fir=Up%Fir
	Up2%Sec=Up%Sec
	Up2%Thir=Up%Thir
	Up2%Four=Up%Four

	do while(abs(error)>10e-7)
	error=0
	Temp%Fir=Up2%Fir
	Temp%Sec=Up2%Sec
	Temp%Thir=Up2%Thir
	Temp%Four=Up2%Four

	Rn1=RU(u,v,p,rho)
    Up2%Fir=Up%Fir+dt*(a21*Rn%Fir+a22*Rn1%Fir)
    Up2%Sec=Up%Sec+dt*(a21*Rn%Sec+a22*Rn1%Sec)
    Up2%Thir=Up%Thir+dt*(a21*Rn%Thir+a22*Rn1%Thir)
    Up2%Four=Up%Four+dt*(a21*Rn%Four+a22*Rn1%Four)

    rho=Jac*Up2%Fir
    u=Jac*Up2%Sec/rho
    v=Jac*Up2%Thir/rho
    p=(gama-1.0)*(Jac*Up2%Four-rho/2.0*(u**2+v**2))
	p(1,:)=p0
	rho(1,:)=rho0
	u(1,:)=u0
	v(1,:)=0

	Up2%Fir=rho/Jac
	Up2%Sec=rho*u/Jac
	Up2%Thir=rho*v/Jac
	Up2%Four=eng/Jac

	do i=1,N1
	do j=1,N2
	 if (abs(Up2%Fir(i,j)-Temp%Fir(i,j))>error) then
	 error=abs(Up2%Fir(i,j)-Temp%Fir(i,j))
	 end if
	if (abs(Up2%Sec(i,j)-Temp%Sec(i,j))>error) then
	 error=abs(Up2%Sec(i,j)-Temp%Sec(i,j))
	 end if
	if (abs(Up2%Thir(i,j)-Temp%Thir(i,j))>error) then
	 error=abs(Up2%Thir(i,j)-Temp%Thir(i,j))
	 end if
	if (abs(Up2%Four(i,j)-Temp%Four(i,j))>error) then
	 error=abs(Up2%Four(i,j)-Temp%Four(i,j))
	 end if
	end do
	end do
	print *,error
	end do

!==================================F3=================================
	error=1
	Up3%Fir=Up%Fir
	Up3%Sec=Up%Sec
	Up3%Thir=Up%Thir
	Up3%Four=Up%Four

	do while(abs(error)>10e-7)
	error=0
	Temp%Fir=Up3%Fir
	Temp%Sec=Up3%Sec
	Temp%Thir=Up3%Thir
	Temp%Four=Up3%Four

	Rn2=RU(u,v,p,rho)
    Up3%Fir=Up%Fir+dt*(a31*Rn%Fir+a32*Rn1%Fir+a33*Rn2%Fir)
    Up3%Sec=Up%Sec+dt*(a31*Rn%Sec+a32*Rn1%Sec+a33*Rn2%Sec)
    Up3%Thir=Up%Thir+dt*(a31*Rn%Thir+a32*Rn1%Thir+a33*Rn2%Thir)
    Up3%Four=Up%Four+dt*(a31*Rn%Four+a32*Rn1%Four+a33*Rn2%Four)

    rho=Jac*Up3%Fir
    u=Jac*Up3%Sec/rho
    v=Jac*Up3%Thir/rho
    p=(gama-1.0)*(Jac*Up3%Four-rho/2.0*(u**2+v**2))
	p(1,:)=p0
	rho(1,:)=rho0
	u(1,:)=u0
	v(1,:)=0

	Up3%Fir=rho/Jac
	Up3%Sec=rho*u/Jac
	Up3%Thir=rho*v/Jac
	Up3%Four=eng/Jac

	do i=1,N1
	do j=1,N2
	 if (abs(Up3%Fir(i,j)-Temp%Fir(i,j))>error) then
	 error=abs(Up3%Fir(i,j)-Temp%Fir(i,j))
	 end if
	if (abs(Up3%Sec(i,j)-Temp%Sec(i,j))>error) then
	 error=abs(Up3%Sec(i,j)-Temp%Sec(i,j))
	 end if
	if (abs(Up3%Thir(i,j)-Temp%Thir(i,j))>error) then
	 error=abs(Up3%Thir(i,j)-Temp%Thir(i,j))
	 end if
	if (abs(Up3%Four(i,j)-Temp%Four(i,j))>error) then
	 error=abs(Up3%Four(i,j)-Temp%Four(i,j))
	 end if
	end do
	end do
	print *,error
	end do

!================================================================



    Up%Fir=Up%Fir+dt*(b1*Rn%Fir+b2*Rn1%Fir+b3*Rn2%Fir)
    Up%Sec=Up%Sec+dt*(b1*Rn%Sec+b2*Rn1%Sec+b3*Rn2%Sec)
    Up%Thir=Up%Thir+dt*(b1*Rn%Thir+b2*Rn1%Thir+b3*Rn2%Thir)
    Up%Four=Up%Four+dt*(b1*Rn%Four+b2*Rn1%Four+b3*Rn2%Four)

    TimePace%Fir=Jac*Up%Fir
    TimePace%Sec=Jac*Up%Sec/TimePace%Fir
    TimePace%Thir=Jac*Up%Thir/TimePace%Fir
    TimePace%Four=(gama-1.0)*(Jac*Up%Four-TimePace%Fir/2.0*(TimePace%Sec**2+TimePace%Thir**2))
	return
end