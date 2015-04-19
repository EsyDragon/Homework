!分子热运动地舍选抽样方法
!f(x)=2/sqrt(PI)*sqrt(x)*exp(-x)
!乘分布f(x)=h(x)f1(x)
!其中f1(x)=2.0/3.0*(exp(-2.0/3.0*x))
!满足概率密度为f1(x)的随机数通过直接取样法实现
!h(x)=(3.0*exp(-x/3.0)*sqrt(x))/sqrt(PI)

      module global
            implicit none
            real,parameter::PI=3.1415926
            real::M
      end module

      program test
            use global
            implicit none
            real::x0=1E-2
            real::r
            integer::i
            
            call random_seed()
            call finM(x0)

            open(100,file='homework0502.dat',status='replace')
            do i=1,10000
                  call thermalMotion(r)
                  write(100,'(F12.8)') r
            end do
            close(100)
      end

      subroutine finM(x0)
            use global
            implicit none
            real::x0
            real,external::h
            real::step=0.01
            integer::i
            real::temp

            M=h(x0)
            do i=1,1000
                  temp=h(x0+i*step)
                  if (M > temp) then
                        exit
                  else
                        M=temp
                  end if
            end do
      end subroutine
                  
      subroutine thermalMotion(x)
            use global
            implicit none
            real::x,r,r1
            real,external::f,h
            do 
                  call F1Distribution(x,r)
                  call random_number(r1)
                  if (h(x)/M >= r1) then
                        exit
                  end if
            end do
      end subroutine

      subroutine F1Distribution(x,r)
            use global
            implicit none
            real::x,r

            call random_number(r)
            x=3.0/2.0*log(1.0/(1.0-r))
      end subroutine

      real function f(x)
            use global
            implicit none
            real::x

            f=2/sqrt(PI)*sqrt(x)*exp(-x)

            return
      end


      real function f1(x)
            use global
            implicit none
            real::x

            f1=2.0/3.0*(exp(-2.0/3.0*x))

            return 
      end

      real function h(x)
            use global
            implicit none
            real::x
            real,external::f
            real,external::f1
            
            h=f(x)/f1(x)

            return
      end

            


