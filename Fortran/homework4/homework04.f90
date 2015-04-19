!单一边界下的边界元法
      include'gauss.F90'
      module global
            implicit none
            integer,parameter::N=32,M=16,Allpoint=N+M
            real::x(N),y(N)
            real::xmid(N),ymid(N)
            real::xinner(M),yinner(M)
            integer,parameter::BoundaryGamma=16
            real::UandQ0(N)
            real::PI=3.1415926
      end module

      program main
            use global
            implicit none
            real::H(Allpoint,N),G(Allpoint,N),B(Allpoint)
            real::UandQ(N),UInner(M)
            integer::i

            !初始化坐标
            do i=1,BoundaryGamma
                  x(i)=10*cos(2.0*PI/BoundaryGamma*(i-1))
                  y(i)=10*sin(2.0*PI/BoundaryGamma*(i-1))
            end do

            do i=(BoundaryGamma+1),N
                  x(i)=6*cos(2.0*PI                               &
     &                        /BoundaryGamma*(i-BoundaryGamma-1))
                  y(i)=6*sin(2.0*PI                               &
     &                        /BoundaryGamma*(i-BoundaryGamma-1))
            end do

            
            do i=1,BoundaryGamma
                  xmid(i)=(x(i)+x(mod(i,BoundaryGamma)+1))/2.0
                  ymid(i)=(y(i)+y(mod(i,BoundaryGamma)+1))/2.0
            end do

            do i=(BoundaryGamma+1),N
                  xmid(i)=(x(i)+x(mod(i,BoundaryGamma)            &
     &                        +BoundaryGamma+1))/2.0
                  ymid(i)=(y(i)+y(mod(i,BoundaryGamma)            &
     &                        +BoundaryGamma+1))/2.0
            end do

            do i=1,M
                  xinner(i)=8*cos(2.0*PI/BoundaryGamma            &
     &                        *(i-BoundaryGamma-1))
                  yinner(i)=8*sin(2.0*PI/BoundaryGamma            &
     &                        *(i-BoundaryGamma-1))

            end do


            !初始化gamma1下的第一类边界和gamma2下的第二类边界
            do i=1,BoundaryGamma
                  UandQ0(i)=100
            end do
            !q
            do i=(BoundaryGamma+1),N
                  UandQ0(i)=-12
            end do
            !计算
            call calculateHGB(H,G,B)
            call calculateBoudary(H,G,B,UandQ)
            call calculateInner(H,G,B,UandQ,UInner)
            !储存数据
            open(100,file='homework04.dat',status='replace')
200         format(F20.8,1X,F20.8,1X,F20.8)
            do i=1,BoundaryGamma
                  write(100,200) xmid(i),ymid(i),UandQ0(i)
            end do
            do i=1,M
                  write(100,200) xinner(i),yinner(i),UInner(i)
            end do
            do i=(BoundaryGamma+1),N
                  write(100,200) xmid(i),ymid(i),UandQ(i)
            end do
            close(100)
      end

      subroutine calculateHGB(H,G,B)
            use global
            implicit none
            integer::i,j
            real::H(Allpoint,N),G(Allpoint,N),B(Allpoint)
            real::s

            do j=1,BoundaryGamma
                  do i=1,N
                        call calculateElement(xmid(i),ymid(i),    &
     &                  x(j),y(j),x(mod(j,BoundaryGamma)+1),      &
     &                  y(mod(j,BoundaryGamma)+1),H(i,j),G(i,j))
                  end do
            end do
            
            do j=(BoundaryGamma+1),N
                  do i=1,N
                        call calculateElement(xmid(i),ymid(i),    &
     &                  x(mod(j,BoundaryGamma)+BoundaryGamma+1),  &
     &                  y(mod(j,BoundaryGamma)+BoundaryGamma+1),  &
     &                  x(j),y(j),                                &
     &                  H(i,j),G(i,j))
                  end do
            end do

            do i=1,N
                  H(i,i)=-PI
                  if (i<=BoundaryGamma) then 
                        s=sqrt((x(mod(i,BoundaryGamma)+1)-x(i))**2&
     &                  +(y(mod(i,BoundaryGamma)+1)-y(i))**2)
                  else
                        s=sqrt((x(mod(i,BoundaryGamma)            &
     &                  +BoundaryGamma+1)-x(i))**2+               &
     &                  (y(mod(i,BoundaryGamma)+BoundaryGamma+1)  &
     &                  -y(i))**2)
                  end if
                  G(i,i)=s*(log(s/2.0)-1.0)
            end do

            do j=1,BoundaryGamma
                  do i=(N+1),Allpoint
                        call calculateElement(xinner(i-N),        &
     &                  yinner(i-N),x(j),y(j),                    &
     &                  x(mod(j,BoundaryGamma)+1),                &
     &                  y(mod(j,BoundaryGamma)+1),H(i,j),G(i,j))
                  end do
            end do
            !inner的H部分
            do j=(BoundaryGamma+1),N
                  do i=(N+1),Allpoint
                        call calculateElement(xinner(i-N),        &
     &                  yinner(i-N),                              &
     &                  x(mod(j,BoundaryGamma)+BoundaryGamma+1),  &
     &                  y(mod(j,BoundaryGamma)+BoundaryGamma+1),  &
     &                  x(j),y(j),                                &
     &                  H(i,j),G(i,j))
                  end do
            end do
            !f=0
            B(:)=0
      end subroutine

      subroutine calculateElement(xi,yi,x1,y1,x2,y2,H,G)
            implicit none
            real::xi,yi
            real::x1,y1,x2,y2
            real::H,G
            real::x21,y21,x1i,y1i,x2i,y2i
            real(kind=8)::s,d,r1,r2,rd,theta

            x21=x2-x1
            y21=y2-y1
            x1i=x1-xi
            y1i=y1-yi
            x2i=x2-xi
            y2i=y2-yi
            
            s=sqrt(x21**2+y21**2)
            d=-(x1i*x21+y1i*y21)/s
            r1=sqrt(x1i**2+y1i**2)
            r2=sqrt(x2i**2+y2i**2)
            rd=(x1i*y21-y1i*x21)/s
            theta=acos((x1i*x2i+y1i*y2i)/(r1*r2))
            H=rd/abs(rd)*theta
            G=(s-d)*log(r2)+d*log(r1)-s+abs(rd)*theta
      end subroutine

      subroutine calculateBoudary(H,G,B,UandQ)
            use global
            use gauss
            implicit none
            real::H(Allpoint,N),G(Allpoint,N),B(Allpoint)
            real::A(N,N),R(N)
            real::UandQ(N)
            integer::i,j

            do i=1,N
                  R(i)=0
            end do
            
            do j=1,N
                  if (j<=BoundaryGamma) then
                        do i=1,N
                              A(i,j)=G(i,j)
                              R(i)=R(i)+H(i,j)*UandQ0(j)
                        end do
                  else
                        do i=1,N
                              A(i,j)=-H(i,j)
                              R(i)=R(i)-G(i,j)*UandQ0(j)
                        end do
                  end if
            end do
            do i=1,N
                  R(i)=R(i)+B(i)
            end do
            call solve(A,R,UandQ,N)
      end subroutine

      subroutine calculateInner(H,G,B,UandQ,UInner)
            use global
            implicit none
            real::H(Allpoint,N),G(Allpoint,N),B(Allpoint)
            real::UInner(M),UandQ(N)
            integer::i,j

            UInner(:)=0

            do j=1,N
                  if (j<=BoundaryGamma) then
                        do i=N+1,Allpoint
                              UInner(i-N)=UInner(i-N)             &
      &                       +H(i,j)*UandQ0(j)-G(i,j)*UandQ(j)
                        end do
                  else
                        do i=N+1,Allpoint
                              UInner(i-N)=UInner(i-N)             &
     &                        +H(i,j)*UandQ(j)-G(i,j)*UandQ0(j)
                        end do
                  end if
            end do

            do i=N+1,Allpoint
                  UInner(i-N)=UInner(i-N)+B(i-N)
            end do
            UInner(:)=UInner(:)/(2*PI)
      end subroutine
