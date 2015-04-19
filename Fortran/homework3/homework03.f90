!变量精度所限COLUMN<82
      module global
            implicit none
            integer,parameter::ROW=3,COLUMN=10
            integer,parameter::E1=(COLUMN-1),E0=2*(COLUMN-1)*(ROW-1)
            integer,parameter::N1=(ROW-1)*COLUMN,N0=ROW*COLUMN
            real,parameter::PI=ACOS(-1.0)
            real::RADIUS(ROW)
            real::X(N0),Y(N0)
            integer::V(E0,3)
            real::Z(E0,3,3),K(N1,N1),R(N0)
            real::rf(E0,3),rq(E0,3)
      end module

      program main
            use global
            implicit none
            integer::i,e,j
            real::u(N0)

            call VNprpr()
            
            do j=1,ROW
                  RADIUS(j)=6.0+(10.0-6.0)/(ROW-1)*(j-1)
            end do

            do j=1,ROW
                  do i=0,E1
                        X(COLUMN*(j-1)+1+i)=RADIUS(j)*cos(PI/2.0/E1*i)
                        Y(COLUMN*(j-1)+1+i)=RADIUS(j)*sin(PI/2.0/E1*i)
                  end do
            end do

            !第一类边界条件
            do i=(N1+1),N0
                  u(i)=100
            end do

            !单元分析
            do e=1,E0
                  call elementAnalysis(e)
            end do
            !计算K,R
            call asemble(u)
            call solve(u)
            open(100,file='homework03.dat',status='replace')
            do i=1,N0
                  write(100,'(F20.8XF20.8XF20.8)') X(i),Y(i),u(i)
            end do
            close(100)
      end

      subroutine VNprpr()
            use global
            implicit none
            integer::e,i,j

            do e=1,E1
                  V(e,1)=COLUMN+e
                  V(e,2)=e+1
                  V(e,3)=e
            end do

            do j=2,2*(ROW-1),2
                  do i=1,(COLUMN-1)
                        e=(COLUMN-1)*(j-1)+i
                        V(e,1)=COLUMN*(j/2-1)+i+1
                        V(e,2)=COLUMN*(j/2)+i
                        V(e,3)=COLUMN*(j/2)+i+1
                  end do
            end do

            do j=3,2*(ROW-1),2
                  do i=1,(COLUMN-1)
                        e=(COLUMN-1)*(j-1)+i
                        V(e,1)=COLUMN*(j/2)+i
                        V(e,2)=COLUMN*(j/2)+i+1
                        V(e,3)=COLUMN*(j/2+1)+i
                  end do
            end do
      end

      subroutine elementAnalysis(e)
            use global
            implicit none
            integer::e
            real::meanF,meanQ
            real::xe(3),ye(3)
            integer::i,j,n
            real,external::f,q
            real::b(3),c(3)
            real::s1,d

            meanF=0
            meanQ=0
            !根据V-n表对x(e,i)=x(n)赋值
            do i=1,3
                  n=V(e,i)
                  xe(i)=X(n)
                  ye(i)=Y(n)
                  meanF=meanF+f(xe(i),ye(i))/3.0
            end do

            if (e <= E1) then
                  meanQ=(q(xe(2),ye(2))+q(xe(3),ye(3)))/2.0
            end if
            
            !计算Rf,Rq,Z
            b(1)=ye(2)-ye(3)
            b(2)=ye(3)-ye(1)
            b(3)=ye(1)-ye(2)
            c(1)=xe(3)-xe(2)
            c(2)=xe(1)-xe(3)
            c(3)=xe(2)-xe(1)

            d=abs(0.5*(b(1)*c(2)-b(2)*c(1)))
                  !Z,Rf
            do i=1,3
                  do j=1,i
                        z(e,i,j)=0.25*(b(i)*b(j)+c(i)*c(j))/d
                        z(e,j,i)=z(e,i,j)
                  end do
                  rf(e,i)=(-1.0/3.0)*d*meanF
            end do
                  !Rq
            rq(e,:)=0
            if (e <= E1) then
                  s1=sqrt(b(1)**2+c(1)**2)
                  rq(e,2)=0.5*s1*meanQ
                  rq(e,3)=0.5*s1*meanQ
            end if
      end subroutine

      subroutine asemble(u)
            use global
            implicit none
            real::u(N0)
            integer::e,i,j,n,m

            !K矩阵的赋值
            K=0
            do e=1,E0
                  do i=1,3
                        n=V(e,i)
                        do j=1,3
                              m=V(e,j)
                              if (n <= N1) then
                                    if (m <= N1) then
                                          K(n,m)=K(n,m)+z(e,i,j)
                                    end if
                              end if
                        end do
                  end do
            end do
            
            !R=R+Rf的赋值
            R=0
            do e=1,E0
                  do i=1,3
                        n=V(e,i)
                        if (n <= N1) then
                              R(n)=R(n)+rf(e,i)
                        end if
                  end do
            end do

            !R=R+Rq的赋值
            do e=1,E1
                  do i=2,3
                        n=V(e,i)
                        if (n <= N1) then
                              R(n)=R(n)+rq(e,i)
                        end if
                  end do
            end do

            !R=R+Ru的赋值
            do e=1,E0
                  do i=1,3
                        n=V(e,i)
                        do j=1,3
                              m=V(e,j)
                              if (m >= (N1+1)) then
                                    R(n)=R(n)-z(e,i,j)*u(m)
                              end if
                        end do
                  end do
            end do
      end subroutine

      subroutine solve(u)
            use global
            implicit none
            real::u(N1)
            real::w,eps=1E-4
            integer::i
            real p,q,temp

            w=2.0-sqrt(2.0)*PI*sqrt(1.0/(N1**2))

            do
                  p=0
                  do i=1,N1
                        temp=u(i)
                        u(i)=W/K(i,i)*(R(i)-sum(K(i,:)*u(:))+K(i,i)*u(i))+(1-W)*u(i)
                        q=abs(u(i)-temp)
                        if (p <= q) p=q
                  end do

                  if(p < eps) then
                        exit
                  end if
            end do

      end subroutine

      real function f(x,y)
            implicit none
            real::x,y

            f=4
            
            return
      end

      real function q(x,y)
            implicit none
            real::x,y

            q=-12

            return
      end


