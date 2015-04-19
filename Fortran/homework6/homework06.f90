!含时薛定谔的解法
      module global
            implicit none
            real(kind=8),parameter::L=1E-4
            real(kind=8),parameter::T=5.17806364251862E-2
            real(kind=8),parameter::Mu=1.673E-27
            real(kind=8),parameter::Omega=200
            real(kind=8),parameter::Hbar=1.055E-34
            real(kind=8),parameter::Alpha=sqrt(Mu*Omega/Hbar)
            real(kind=8),parameter::PI=3.1415926535897
            complex(kind=8),parameter::complex_I=(0,1)
      end module

      program main
            use global
            implicit none
            real(kind=8),external::U
            real(kind=8),external::Psi
            real(kind=8),external::Integrate
            integer,parameter::Kmax=10
            integer,parameter::Xmax=40
            integer,parameter::Tmax=40
            real(kind=8)::deltaX=2*L/(Xmax-1),deltaT=T/(Tmax-1)
            real(kind=8)::c(Kmax),E(Kmax)
            complex(kind=8)::temp(Xmax,Tmax)
            real(kind=8)::squarePsi(Xmax,Tmax)
            integer::i,j,k,n

            temp(:,:)=0.0

            do k=1,Kmax
                  n=2*k-1
                  c(k)=Integrate(U,n,Psi,-L,L)
                  E(k)=((PI*Hbar*n)**2)/(8*Mu*L**2)
            end do


            do j=1,Tmax
                  do i=1,Xmax
                        do k=1,Kmax
                              temp(i,j)=temp(i,j)+c(k)*exp(-complex_I*E(k)*(deltaT*(j-1))/Hbar)*U(2*k-1,-L+deltaX*(i-1))
                        end do
                        squarePsi(i,j)=abs(temp(i,j))**2
                  end do
            end do

            open(100,file="homework06.dat")
            do j=1,Tmax
                  do i=1,Xmax
                        write(100,'(F20.8XF20.8XF20.8)') -L+deltaX*(i-1),deltaT*(j-1),squarePsi(i,j)
                  end do
            end do
            close(100)
            end

      real(kind=8) function U(n,x)
            use global
            implicit none
            integer::n
            real(kind=8)::x

            U=sin((n*PI*(x+L))/(2.0*L))/sqrt(L)

            return
      end

      real(kind=8) function Psi(x)
            use global
            implicit none
            real(kind=8)::x

            Psi=sqrt(Alpha/sqrt(PI))*exp(-(Alpha*x)**2/2.0)

            return
      end

      real(kind=8) function Integrate(U,n,Psi,a,b)
            implicit none
            real(kind=8)::a,b
            real(kind=8),external::U
            real(kind=8),external::Psi
            integer::n
            integer,parameter::M=100
            real(kind=8)::d
            integer::i
            


            d=(b-a)/M
            Integrate=1.0/2.0*(U(n,a)*Psi(a)+U(n,b)*Psi(b))*d
            do i=1,(M-1)
                  Integrate=Integrate+U(n,a+i*d)*Psi(a+i*d)*d
            end do

            return
      end

            
            
