      subroutine vcalir(ndim, gamma, ifirst, ilast, udemi, cdemi, ecdemi, &
                        a, Lip12 )

      IMPLICIT NONE

      INTEGER :: i, m
      INTEGER :: ndim, ifirst, ilast

      REAL(kind=8) :: gamma, guvw2, usnc
      REAL(kind=8), DIMENSION(3) :: a
      REAL(kind=8), DIMENSION(ifirst:ilast) :: cdemi, ecdemi
      REAL(kind=8), DIMENSION(3, ifirst:ilast) :: udemi

      REAL(kind=8), DIMENSION(ifirst:ilast,5,5) :: Lip12

!
      if(ndim == 1) then

        do i=ifirst,ilast

          guvw2=(gamma-1.)*ecdemi(i)

          usnc = udemi(1,i)*cdemi(i)

          Lip12(i,1,1)=cdemi(i)**2-guvw2
          Lip12(i,2,1)=.5*(guvw2-usnc)
          Lip12(i,3,1)=.5*(guvw2+usnc)

          Lip12(i,1,2)=(gamma-1.)*udemi(1,i)
          Lip12(i,2,2)=-.5*((gamma-1.)*udemi(1,i)-cdemi(i))
          Lip12(i,3,2)=-.5*((gamma-1.)*udemi(1,i)+cdemi(i))

          Lip12(i,1,3) = -(gamma-1.)
          Lip12(i,2,3) = (gamma-1.)*.5
          Lip12(i,3,3) = (gamma-1.)*.5

        enddo

      elseif(ndim == 2) then

        do i=ifirst,ilast

          guvw2=(gamma-1.)*ecdemi(i)

          usnc = 0.
          do m = 1, ndim
            usnc = usnc + a(m)*udemi(m,i)*cdemi(i)
          enddo

          Lip12(i,1,1)=cdemi(i)**2-guvw2
          Lip12(i,2,1)=(a(2)*udemi(1,i) - a(1)*udemi(2,i))*cdemi(i)
          Lip12(i,3,1)=.5*(guvw2-usnc)
          Lip12(i,4,1)=.5*(guvw2+usnc)

          Lip12(i,1,2)=(gamma-1.)*udemi(1,i)
          Lip12(i,2,2)=-a(2)*cdemi(i)
          Lip12(i,3,2)=-.5*((gamma-1.)*udemi(1,i)-a(1)*cdemi(i))
          Lip12(i,4,2)=-.5*((gamma-1.)*udemi(1,i)+a(1)*cdemi(i))

          Lip12(i,1,3)=(gamma-1.)*udemi(2,i)
          Lip12(i,2,3)= a(1)*cdemi(i)
          Lip12(i,3,3)=-.5*((gamma-1.)*udemi(2,i)-a(2)*cdemi(i))
          Lip12(i,4,3)=-.5*((gamma-1.)*udemi(2,i)+a(2)*cdemi(i))

          Lip12(i,1,4) = -(gamma-1.)
          Lip12(i,2,4) = 0.
          Lip12(i,3,4) = (gamma-1.)*.5
          Lip12(i,4,4) = (gamma-1.)*.5

        enddo

      else

        do i=ifirst,ilast

          guvw2=(gamma-1.)*ecdemi(i)

          usnc = 0.
          do m = 1, ndim
            usnc = usnc + a(m)*udemi(m,i)*cdemi(i)
          enddo

          Lip12(i,1,1)=cdemi(i)**2-guvw2
          Lip12(i,2,1)=((a(3)-a(2))*udemi(1,i) + &
                        (a(1)-a(3))*udemi(2,i) + &
                        (a(2)-a(1))*udemi(3,i))*cdemi(i)
          Lip12(i,3,1)=((a(3)**2-a(1)*a(2))*udemi(1,i) + &
                        (a(1)**2-a(2)*a(3))*udemi(2,i) + &
                        (a(2)**2-a(1)*a(3))*udemi(3,i))*cdemi(i)
          Lip12(i,4,1)=.5*(guvw2-usnc)
          Lip12(i,5,1)=.5*(guvw2+usnc)

          Lip12(i,1,2)=(gamma-1.)*udemi(1,i)
          Lip12(i,2,2)=(a(2)-a(3))*cdemi(i)
          Lip12(i,3,2)=(a(1)*a(2)-a(3)**2)*cdemi(i)
          Lip12(i,4,2)=-.5*((gamma-1.)*udemi(1,i)-a(1)*cdemi(i))
          Lip12(i,5,2)=-.5*((gamma-1.)*udemi(1,i)+a(1)*cdemi(i))

          Lip12(i,1,3)=(gamma-1.)*udemi(2,i)
          Lip12(i,2,3)=(a(3)-a(1))*cdemi(i)
          Lip12(i,3,3)=(a(2)*a(3)-a(1)**2)*cdemi(i)
          Lip12(i,4,3)=-.5*((gamma-1.)*udemi(2,i)-a(2)*cdemi(i))
          Lip12(i,5,3)=-.5*((gamma-1.)*udemi(2,i)+a(2)*cdemi(i))

          Lip12(i,1,4)=(gamma-1.)*udemi(3,i)
          Lip12(i,2,4)=(a(1)-a(2))*cdemi(i)
          Lip12(i,3,4)=(a(1)*a(3)-a(2)**2)*cdemi(i)
          Lip12(i,4,4)=-.5*((gamma-1.)*udemi(3,i)-a(3)*cdemi(i))
          Lip12(i,5,4)=-.5*((gamma-1.)*udemi(3,i)+a(3)*cdemi(i))

          Lip12(i,1,5) = -(gamma-1.)
          Lip12(i,2,5) = 0.
          Lip12(i,3,5) = 0.
          Lip12(i,4,5) = (gamma-1.)*.5
          Lip12(i,5,5) = (gamma-1.)*.5

        enddo

      endif

      end subroutine vcalir
