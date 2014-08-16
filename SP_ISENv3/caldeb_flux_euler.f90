      subroutine caldeb_flux_euler(maille, dir, fmp_gauche,f_gauche)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      TYPE(structure_maille), pointer :: maille

      INTEGER :: dir

      REAL(kind=8), DIMENSION(nvar) :: fmp_gauche,f_gauche

      INTEGER :: l, m

      REAL(kind=8), DIMENSION(nvar, 3) :: u_n

      REAL(kind=8) :: rho_ec
      REAL(kind=8), DIMENSION(3) :: psgM2, tsgM2

      INTEGER, DIMENSION(4,3) :: ijk_dt

      REAL(kind=8) :: rho_wall, t_local, p_local
      REAL(kind=8) :: divrhoV
      REAL(kind=8), DIMENSION(ndim, ndim) :: rhoV_xyz

      REAL(kind=8), DIMENSION(1) :: cp12, hp12, ecp12
      REAL(kind=8), DIMENSION(3, 1) :: up12
      REAL(kind=8), DIMENSION(nvar) :: vp

      TYPE(structure_maille), pointer :: courant

! ======== Flux Euler : Condition a la limite amont du domaine ========

! +++++condition de frontiere libre

      if( caltype_deb(dir) == 0 ) then

! ! *****Probleme de Riemann
!
! ! .....Moyenne sur la frontiere
!
!         do m = 2, nvar-1
!           up12(m-1,1) = 0.5*(maille%u(m)/ maille%u(1) + v_left(m-1))
!         enddo
!
! ! .....valeurs propres
!         do m = 1, nvar
!           vp(m) = up12(dir,1)
!         enddo
!
! ! *****Projection des flux Euler sur les vecteurs propres
!
!          fmp_gauche(:) = 0.
!          f_gauche(:) = 0.
!
!           if( vp(1) < 0. ) then
!             fmp_gauche(:) = maille%fmp_dt(:)
!             f_gauche(:) = maille%f_dt(:)
!           endif

! *****Flux Euler nuls a la paroi
        do m = 1, nvar
          fmp_gauche(m) = 0.
          f_gauche(m) = 0.
        enddo

! +++++condition de paroi

      elseif( caltype_deb(dir) == 1 ) then

! *****Flux Euler nuls a la paroi
        do m = 1, nvar
          fmp_gauche(m) = 0.
          f_gauche(m) = 0.
        enddo

! ! *****Extrapolation de la pression (ordre 2)
!
! ! .....1er point au dessus de la paroi
!         courant => maille
!
!         u_n(1,1) = courant%u(1)
!         rho_ec = 0.
!         do m = 2,nvar-1
!           u_n(m,1) = courant%u(m)
!           rho_ec = rho_ec + courant%u(m)**2
!         enddo
!         rho_ec = 0.5 * rho_ec / courant%u(1)
!         u_n(nvar,1) = courant%u(nvar)
!
!         psgM2(1) = (gamma-1.)*(courant%u(nvar)-rho_ec)
!
!         tsgM2(1) = psgM2(1) / courant%u(1)
!
! ! .....2eme point au dessus de la paroi
!         courant => courant%support(dir)%suivant
!
!         u_n(1,2) = courant%u(1)
!         rho_ec = 0.
!         do m = 2,nvar-1
!           u_n(m,2) = courant%u(m)
!           rho_ec = rho_ec + courant%u(m)**2
!         enddo
!         rho_ec = 0.5 * rho_ec / courant%u(1)
!         u_n(nvar,2) = courant%u(nvar)
!
!         psgM2(2) = (gamma-1.)*(courant%u(nvar)-rho_ec)
!
!         tsgM2(2) = psgM2(2) / courant%u(1)
!
! ! .....3eme point au dessus de la paroi
!         courant => courant%support(dir)%suivant
!
!         u_n(1,3) = courant%u(1)
!         rho_ec = 0.
!         do m = 2,nvar-1
!           u_n(m,3) = courant%u(m)
!           rho_ec = rho_ec + courant%u(m)**2
!         enddo
!         rho_ec = 0.5 * rho_ec / courant%u(1)
!         u_n(nvar,3) = courant%u(nvar)
!
!         psgM2(3) = (gamma-1.)*(courant%u(nvar)-rho_ec)
!
!         tsgM2(3) = psgM2(3) / courant%u(1)
!
! ! *****Gradients de vitesses
!         do l = 1, ndim
!           do m = 2, nvar-1
!             rhoV_xyz(m-1, l) = 0.
!           enddo
!         enddo
!
!         do m = 2, nvar-1
!           rhoV_xyz(m-1, dir) = (225.*u_n(m, 1) - 50.*u_n(m, 2) &
!                                + 9.*u_n(m,3))/(60.*maille%dx(dir))
! !!          rhoV_xyz(m-1, dir) = (9.*u_n(m, 1) - u_n(m, 2))/(3.*maille%dx(dir))
!         enddo
!
! ! *****Divergence de vitesse
!         divrhoV = 0.
!         do m = 2, nvar-1
!           divrhoV = divrhoV + rhoV_xyz(m-1,m-1)
!         enddo
!
! ! *****Temperature de paroi
!
!         if( cle_twall ) then
! ! .....paroi adiabatique ==> d/dn Temperature = 0 - ordre 3 (ou ordre 2)
!           t_local = ( 225.*tsgM2(1) - 50.*tsgM2(2) + 9.*tsgM2(3)) / 184.
! !!          t_local = ( 9.*tsgM2(1) - tsgM2(2) ) / 8.
!         else
! ! .....paroi thermostatee
!           t_local = t_wall/(gamma*Mach**2)
!         endif
!
! ! *****Densite a la paroi par equation de continuite
!         rho_wall = ( 15.*u_n(1,1) - 10.*u_n(1,2) + 3.*u_n(1,3) )/8. - &
!                    dt * divrhoV
!
! ! *****Flux Euler suivant (dir) :
! !      pression paroi extrapolee - ordre 3
!
!         p_local = ( 15.*psgM2(1) - 10.*psgM2(2) + 3.*psgM2(3) )/8.
!         fmp_gauche(dir+1) = p_local
!!
!!
!!!      Pression paroi ==> d/dn pression = 0 - ordre 3
!!        fmp_gauche(dir+1) = ( 225.*psgM2(1) - 50.*psgM2(2) + 9.*psgM2(3)) &
!!                            / 184.
!!!      pression = rho * T
!!        fmp_gauche(dir+1) = rho_wall * t_local
         !
      endif

      end subroutine caldeb_flux_euler