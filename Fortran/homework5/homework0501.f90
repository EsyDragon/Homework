!蒲丰投针计算Pi
      program main
            real::d=4,l=3
            real::y,theta
            real(kind=8)::pi=0
            integer::time=0
            integer,parameter::N=1000

            call random_seed() 
            do i=1,N
                  call random_number(y)
                  y=d*y
                  call random_angle(theta)
                  if (y < l*sin(theta)) then
                        time=time+1
                  end if
            end do

            pi=(2*l)/(real(time)/N*d)
            write(*,*) pi
      end

      subroutine random_angle(theta)
            real::theta
            real::x,y
            do 
                  call random_number(x)
                  call random_number(y)
                  if ((x**2+y**2) < 1) then
                        exit
                  end if
            end do
            theta=2*atan2(y,x)
      end subroutine


