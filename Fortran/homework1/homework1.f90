!此程序采用兼用源码格式,在gfortran下测试通过
!张亦龙 20150306
!三点电荷电场强度的计算
!由于fortran的绘图库大多不够强,本程序导出数据并在Mathematica下绘图
      module global
            implicit none
            integer,parameter::size=256
      end module

      program main
            use global
            implicit none
            real::field(size,size)
            integer::i,j

            call calculateField(field)

            open(100,file="homework1.dat",status="replace")
            do j=1,size
                  write(100,"(256(E12.7X))") field(:,j)
            end do

      end

      subroutine calculateField(field)
            use global
            implicit none
            real::field(size,size)
            integer::i,j
            real::x(size),y(size)

            do i=1,size
                  !防止在点电荷处的无穷大,对网格进行了小量偏移
                  x(i)=-2.0+4.0/size*i+0.001
                  y(i)=-2.0+4.0/size*i+0.001
            end do

            do j=1,size
                  do i=1,size
                        field(i,j)= Sqrt((-(x(i)/(x(i)**2                &
     &                  + (-1 + y(j))**2)**1.5) - x(i)/(x(i)**2          &
     &                  + y(j)**2)**1.5 + (2*x(i))/(x(i)**2              &
     &                  + (1 + y(j))**2)**1.5)**2                        &
     &                  +  (-((-1 + y(j))/(x(i)**2                       &
     &                  + (-1 + y(j))**2)**1.5)                          & 
     &                  - y(j)/(x(i)**2 + y(j)**2)**1.5                  &
     &                  + (2*(1 + y(j)))/(x(i)**2                        &
     &                  + (1 + y(j))**2)**1.5)**2)
                  end do
            end do
      end


