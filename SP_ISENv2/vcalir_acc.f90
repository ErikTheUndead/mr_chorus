      subroutine vcalir_acc(ifirst, ilast,rho, cdemi,Lip12, Rip12 )

      IMPLICIT NONE

      INTEGER :: i
      INTEGER :: ifirst, ilast

      REAL(kind=8), DIMENSION(ifirst:ilast) :: cdemi, rho

      REAL(kind=8), DIMENSION(ifirst:ilast,2,2) :: Lip12, Rip12

        do i=ifirst,ilast

				Rip12(i,1,1) = 1.
				Rip12(i,1,2) = 1.
				Rip12(i,2,1) = + cdemi(i)/rho(i)
				Rip12(i,2,2) = - cdemi(i)/rho(i)

				Lip12(i,1,1) = 0.5
				Lip12(i,1,2) = + 0.5*rho(i)/cdemi(i)
				Lip12(i,2,1) = 0.5
				Lip12(i,2,2) = - 0.5*rho(i)/cdemi(i)

        enddo

      end subroutine vcalir_acc
