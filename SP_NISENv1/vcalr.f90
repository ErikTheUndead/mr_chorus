      subroutine vcalr(ndim, ifirst, ilast, udemi, cdemi, ecdemi, &
                       hdemi, a, Rip12 )

      IMPLICIT NONE

      INTEGER :: i, m
      INTEGER :: ndim, ifirst, ilast

      REAL(kind=8) :: uc2, uac, usnc
      REAL(kind=8) :: r22, r32, r42, r23, r33, r43
      REAL(kind=8), DIMENSION(3) :: a
      REAL(kind=8), DIMENSION(ifirst:ilast) :: cdemi, ecdemi, hdemi
      REAL(kind=8), DIMENSION(3, ifirst:ilast) :: udemi

      REAL(kind=8), DIMENSION(ifirst:ilast,5,5) :: Rip12

!
      if(ndim == 1) then

        do i=ifirst,ilast

          uc2=1./cdemi(i)**2

          uac=-1./cdemi(i)

          usnc = udemi(1,i)*cdemi(i)

          Rip12(i,1,1)=uc2
          Rip12(i,2,1)=udemi(1,i)*uc2
          Rip12(i,3,1)=ecdemi(i)*uc2

          Rip12(i,1,2)=uc2
          Rip12(i,2,2)=(udemi(1,i)+cdemi(i))*uc2
          Rip12(i,3,2)=(hdemi(i)+usnc)*uc2

          Rip12(i,1,3)=uc2
          Rip12(i,2,3)=(udemi(1,i)-cdemi(i))*uc2
          Rip12(i,3,3)=(hdemi(i)-usnc)*uc2

        enddo


      elseif(ndim == 2) then

        do i=ifirst,ilast

          uc2=1./cdemi(i)**2

          uac=1./((a(1)*a(2)-1.)*cdemi(i))

          usnc = 0.
          do m = 1, ndim
            usnc = usnc + a(m)*udemi(m,i)*cdemi(i)
          enddo

          Rip12(i,1,1)=uc2
          Rip12(i,2,1)=udemi(1,i)*uc2
          Rip12(i,3,1)=udemi(2,i)*uc2
          Rip12(i,4,1)=ecdemi(i)*uc2

          Rip12(i,1,2)=0.
          Rip12(i,2,2)=-a(2)/cdemi(i)
          Rip12(i,3,2)= a(1)/cdemi(i)
          Rip12(i,4,2)=Rip12(i,2,2)*udemi(1,i)+Rip12(i,3,2)*udemi(2,i)

          Rip12(i,1,3)=uc2
          Rip12(i,2,3)=(udemi(1,i)+a(1)*cdemi(i))*uc2
          Rip12(i,3,3)=(udemi(2,i)+a(2)*cdemi(i))*uc2
          Rip12(i,4,3)=(hdemi(i)+usnc)*uc2

          Rip12(i,1,4)=uc2
          Rip12(i,2,4)=(udemi(1,i)-a(1)*cdemi(i))*uc2
          Rip12(i,3,4)=(udemi(2,i)-a(2)*cdemi(i))*uc2
          Rip12(i,4,4)=(hdemi(i)-usnc)*uc2

        enddo

      else

        do i=ifirst,ilast

          uc2=1./cdemi(i)**2

          uac=1./((a(1)*a(2)+a(1)*a(3)+a(2)*a(3)-1.)*cdemi(i))

          usnc = 0.
          do m = 1, ndim
            usnc = usnc + a(m)*udemi(m,i)*cdemi(i)
          enddo

          Rip12(i,1,1)=uc2
          Rip12(i,2,1)=udemi(1,i)*uc2
          Rip12(i,3,1)=udemi(2,i)*uc2
          Rip12(i,4,1)=udemi(3,i)*uc2
          Rip12(i,5,1)=ecdemi(i)*uc2

          r22=(a(2)*(a(1)*a(3)-a(2)*a(2))-a(3)*(a(2)*a(3)-a(1)*a(1)))*uac
          r32=(a(3)*(a(1)*a(2)-a(3)*a(3))-a(1)*(a(1)*a(3)-a(2)*a(2)))*uac
          r42=(a(1)*(a(2)*a(3)-a(1)*a(1))-a(2)*(a(1)*a(2)-a(3)*a(3)))*uac

          Rip12(i,1,2)=0.
          Rip12(i,2,2)=r22
          Rip12(i,3,2)=r32
          Rip12(i,4,2)=r42
          Rip12(i,5,2)=r22*udemi(1,i)+r32*udemi(2,i)+r42*udemi(3,i)

          r23=(a(3)*a(3)-a(1)*a(3)-a(1)*a(2)+a(2)*a(2))*uac
          r33=(a(1)*a(1)-a(1)*a(2)-a(2)*a(3)+a(3)*a(3))*uac
          r43=(a(2)*a(2)-a(2)*a(3)-a(1)*a(3)+a(1)*a(1))*uac

          Rip12(i,1,3)=0.
          Rip12(i,2,3)=r23
          Rip12(i,3,3)=r33
          Rip12(i,4,3)=r43
          Rip12(i,5,3)=r23*udemi(1,i)+r33*udemi(2,i)+r43*udemi(3,i)

          Rip12(i,1,4)=uc2
          Rip12(i,2,4)=(udemi(1,i)+a(1)*cdemi(i))*uc2
          Rip12(i,3,4)=(udemi(2,i)+a(2)*cdemi(i))*uc2
          Rip12(i,4,4)=(udemi(3,i)+a(3)*cdemi(i))*uc2
          Rip12(i,5,4)=(hdemi(i)+usnc)*uc2

          Rip12(i,1,5)=uc2
          Rip12(i,2,5)=(udemi(1,i)-a(1)*cdemi(i))*uc2
          Rip12(i,3,5)=(udemi(2,i)-a(2)*cdemi(i))*uc2
          Rip12(i,4,5)=(udemi(3,i)-a(3)*cdemi(i))*uc2
          Rip12(i,5,5)=(hdemi(i)-usnc)*uc2

        enddo

      endif

      end subroutine vcalr
