      !有限元差分法求解矩形域的泊松方程
      !计算物理学p182
      !张亦龙 12344045
      module global
            implicit none
            !W超松弛因子,EPS误差界限,H步长
            real,parameter::W=1.25,EPS=1E-4,H=0.2,PI=3.1415926
            integer,parameter::m=2/H+1,n=1/H+1
      end module
      
      program main
            use global
            implicit none
            integer::i,j
            real::u(m,n),t,p,q
            real::x(m),y(n)
            integer::count=1
            !边界条件
            do i=1,m
                  x(i)=(i-1)*H
                  u(i,1)=x(i)
                  u(i,n)=x(i)*exp(1.0)
            end do
            do j=1,n
                  y(j)=(j-1)*h
                  u(1,j)=0
                  u(m,j)=2*exp(y(j))
            end do
            !初值
            do j=2,n-1
                  do i=2,m-1
                        u(i,j)=0.0
                  end do
            end do
            do
                  !p为每次迭代的最大误差
                  p=0
                  do j=2,n-1
                        do i=2,m-1
                              count=count+1
                              t=u(i,j)
                              u(i,j)=W*(u(i-1,j)+u(i+1,j)+u(i,j-1)             &
     &                        +u(i,j+1)-H**2*x(i)*exp(y(j)))/4                 &
     &                        +(1-W)*u(i,j)
                              q=abs(u(i,j)-t)
                              if (p<=q) p=q
                        end do
                  end do
                  if (p<eps) then
                        exit
                  end if
            end do
            !注意gfortran中当recl=n时,为每个字段为n byte!对于real为4
            open(100,file='homework0201.dat',status='replace',form='unformatted',       &
     &      access='direct',recl=4*m)
            do j=1,n
                  write(100,rec=j) u(:,j)
            end do
            close(100)
            write(*,*) count
      end
      
