      subroutine flux_euler(dir, ifirst, ilast, u_n, feuler, gh, temperature )

      use mod_common

      IMPLICIT NONE

! +++++Variables arguments

      INTEGER, INTENT(IN) :: dir
      INTEGER, INTENT(IN) :: ifirst, ilast

      REAL(kind=8), INTENT(IN), DIMENSION(nvar, ifirst:ilast) :: u_n
      REAL(kind=8), INTENT(OUT), DIMENSION(nvar, ifirst:ilast) ::  feuler
      REAL(kind=8), INTENT(OUT), DIMENSION(ifirst:ilast) :: gh, temperature

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
        rho_ec = 0.5 * rho_ec / u_n(1,l)

! .....Pression
        pression = gamma*(gamma-1.) * ( u_n(nvar,l) - rho_ec ) * Mach**2

! .....Enthalpie totale
        gh(l) = ( u_n(nvar,l) + pression/(gamma*Mach**2) ) / u_n(1,l)

! .....Temperature statique
        temperature(l) = pression / u_n(1,l)

! +++++Flux Euler

        feuler(1,l) = u_n(dir+1,l)

        do m = 2, nvar-1
          feuler(m,l) = u_n(m,l)*feuler(1,l)/u_n(1,l)
        enddo
        feuler(dir+1,l) =  feuler(dir+1,l) + pression/(gamma*Mach**2)

        feuler(nvar,l) = feuler(1,l)*gh(l)

      enddo

      end subroutine flux_euler
