!p173页例题
!一维扩散方程求解
      module global
            implicit none
            real,parameter::H=0.1,TAO=1E-1,D=0.1,TMAX=10
            real,parameter::A1=1,B1=1,C1=0,XMIN=0
            real,parameter::A2=1,B2=-1,C2=0,XMAX=1.0
            integer,parameter::N=1+(XMAX-XMIN)/H
            integer,parameter::KMAX=TMAX/TAO
            real,parameter::P=TAO*D/(H**2)
            real::P1=1/P+1,P2=1/P-1
            integer::s
            real::X(N)=(/ ((XMIN+(s-1)*H),s=1,N) /)
      end module

      program main
            use global
            implicit none
            ! a*u=r
            real::a(N,3)
            real::r(N)
            integer::i,k
            real::t
            !初始条件
            real::u(N)=(/ ((exp(x(I))),I=1,N) /)

            !计算原始系数矩阵
            call coef(a)
            !计算消元系数矩阵
            call eliminationCoef(a)
            
            open(100,file='homework0202.dat',status='replace')
            do k=1,KMAX
                  t=k*TAO
                  
                  r(1)=(B1*P2-H*A1)*u(1)+B1*u(2)+2*H*C1
                  r(N)=B2*u(N-1)+(B2*P2-H*A2)*u(N)+2*H*C2
                  do i=2,N-1
                        r(i)=u(i-1)+2*P2*u(i)+u(i+1)
                  end do
                  call solve(a,u,r)
                  !写入文件
                  do i=1,N
                        write(100,*) x(i),t,u(i)
                  end do
            end do
            close(100)
      end

      subroutine coef(a)
            use global
            implicit none
            real::a(N,3)
            integer::i

            a(1,2)=P1*B1+H*A1
            a(1,3)=-B1
            a(N,1)=-B2
            a(N,2)=B2*P1+H*A2
            do i=2,N-1
                  a(i,1)=-1
                  a(i,2)=2*P1
                  a(i,3)=-1
            end do
      end

      subroutine eliminationCoef(a)
            use global
            implicit none
            integer::i
            real::a(N,3)
            do i=2,N
                  a(i,2)=a(i,2)-a(i,1)*a(i-1,3)/a(i-1,2)
            end do
      end

      subroutine solve(a,u,r)
            use global
            implicit none
            real::a(N,3)
            real::u(N)
            real::r(N)
            integer::i
            do i=2,N
                  r(i)=r(i)-a(i,1)*r(i-1)/a(i-1,2)
            end do
            
            u(N)=r(N)/a(n,2)
            do i=N-1,1,-1
                  u(i)=(r(i)-a(i,3)*u(i+1))/a(i,2)
            end do
      end
