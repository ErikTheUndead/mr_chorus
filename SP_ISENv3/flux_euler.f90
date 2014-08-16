      subroutine flux_euler(ifirst, ilast, u_n, gh)

      use mod_common

      IMPLICIT NONE

! +++++Variables arguments

      INTEGER, INTENT(IN) :: ifirst, ilast

      REAL(kind=8), INTENT(IN), DIMENSION(nvar, ifirst:ilast) :: u_n
      REAL(kind=8), INTENT(OUT), DIMENSION(ifirst:ilast) :: gh

! +++++Variables locales

      INTEGER :: l, m

      REAL(kind=8) :: rho_ec, pression

! ======================= Flux Euler =============================

      do l = ifirst, ilast

! .....Energie cinetique
        rho_ec = 0.
        do m = 2,nvar-1
          rho_ec = rho_ec + u_n(m,l)**2
        enddo
        rho_ec = 0.5 * rho_ec * u_n(1,l)

! .....Pression
        pression = gamma*(gamma-1.) * ( u_n(nvar,l)*u_n(1,l) - rho_ec ) * Mach**2

! .....Enthalpie totale
        gh(l) = ( u_n(nvar,l)*u_n(1,l) + pression/(gamma*Mach**2) ) / u_n(1,l)

      end do

      end subroutine flux_euler
