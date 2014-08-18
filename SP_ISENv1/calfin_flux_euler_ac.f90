      subroutine calfin_flux_euler_ac(maille, dir, fmp_droit,f_droit)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      TYPE(structure_maille), pointer :: maille

      INTEGER :: dir

      REAL(kind=8), DIMENSION(2) :: fmp_droit,f_droit
      REAL(kind=8), DIMENSION(2) :: df,dfmp

      INTEGER :: l, m

      REAL(kind=8), DIMENSION(nvar, 3) :: u_n

      REAL(kind=8) :: rho_ec
      REAL(kind=8), DIMENSION(3) :: psgM2, tsgM2

      INTEGER, DIMENSION(4,3) :: ijk_dt

      REAL(kind=8) :: rho_wall, t_local, p_local
      REAL(kind=8) :: divrhoV
      REAL(kind=8), DIMENSION(ndim, ndim) :: rhoV_xyz


      REAL(kind=8), DIMENSION(1) :: cp12, hp12, ecp12, rho12
      REAL(kind=8), DIMENSION(3, 1) :: up12
      REAL(kind=8), DIMENSION(2) :: vp
      REAL(kind=8), DIMENSION(1,2,2) :: Rip12, Lip12

      TYPE(structure_maille), pointer :: courant

! ======== Flux Euler : Condition a la limite aval du domaine ========

!       print*,'calfin_flux_ac_debut droite', f_droit(:), fmp_droit(:)

! +++++condition de frontiere libre

      if( caltype_fin(dir) == 0 ) then

! *****Probleme de Riemann
! .....Moyenne de Roe sur la frontiere

        ecp12(1) = 0.
        do m = 2, nvar-1
          up12(m-1,1) = 0.5*(maille%u(m)/ maille%u(1) + v_right(m-1))
          ecp12(1) = ecp12(1) + up12(m-1,1)**2
        enddo
        ecp12(1) = 0.5 * ecp12(1)

        hp12(1) = gamma * maille%u(nvar)/maille%u(1) - &
                  (gamma-1.)*ecp12(1)

        cp12(1) = sqrt( (gamma-1.)*(hp12(1) - ecp12(1)) )

        rho12 = 0.5*(maille%u(1)+rho_right)

! .....vecteurs propres : matrice a gauche
        call vcalir_acc(1, 1, rho12, cp12, Lip12, Rip12)

! .....valeurs propres

        vp(1) = + cp12(1)
        vp(2) = - cp12(1)

! *****Projection des flux Euler sur les vecteurs propres

        do m = 1, 2
			df(m) = 0.
         dfmp(m) = 0.
			if( vp(m) > 0. ) then
            do l = 1, 2
              df(m) = df(m) + Lip12(1,m,l)*( maille%f_dt(l) )
!               dfmp(m) = dfmp(m) + Lip12(1,m,l)*( maille%fmp_dt(l) )
            enddo
          endif
        enddo

! *****Projection des flux Euler dans l'espace physique

        do m = 1, 2
          f_droit(m) = 0.
          fmp_droit(m) = 0.
          do l = 1, 2
            f_droit(m) = f_droit(m) + Rip12(1,m,l) * df(l)
!             fmp_droit(m) = fmp_droit(m) + Rip12(1,m,l) * dfmp(l)
          enddo
        enddo

! +++++condition de paroi

      elseif( caltype_fin(dir) == 1 ) then

! *****Probleme de Riemann
! .....Moyenne de Roe sur la frontiere

        ecp12(1) = 0.
        do m = 2, nvar-1
          up12(m-1,1) = v_left(m-1)
          ecp12(1) = ecp12(1) + up12(m-1,1)**2
        enddo
        ecp12(1) = 0.5 * ecp12(1)

        hp12(1) = gamma * maille%u(nvar)/maille%u(1) - &
                  (gamma-1.)*ecp12(1)

        cp12(1) = sqrt( (gamma-1.)*(hp12(1) - ecp12(1)) )

        rho12 = maille%u(1)

! .....vecteurs propres : matrice a gauche
        call vcalir_acc(1, 1, rho12, cp12, Lip12, Rip12)

! .....valeurs propres

        vp(1) = + cp12(1)
        vp(2) = - cp12(1)

! *****Projection des flux Euler sur les vecteurs propres

        do m = 1, 2
			df(m) = 0.
			if( vp(m) > 0. ) then
            do l = 1, 2
              df(m) = df(m) + Lip12(1,m,l)*( maille%f_dt(l) )
!               dfmp(m) = dfmp(m) + Lip12(1,m,l)*( maille%fmp_dt(l) )
            enddo
          endif
        enddo

! *****Projection des flux Euler dans l'espace physique

        do m = 1, 2
          f_droit(m) = 0.
          fmp_droit(m) = 0.
          do l = 1, 2
            f_droit(m) = f_droit(m) + Rip12(1,m,l) * df(l)
!             fmp_droit(m) = fmp_droit(m) + Rip12(1,m,l) * dfmp(l)
          enddo
        enddo


      endif

!       print*,'calfin_flux_ac_fin droite', f_droit(:), fmp_droit(:)

      end subroutine calfin_flux_euler_ac
