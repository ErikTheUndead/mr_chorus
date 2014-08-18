      subroutine residu

      use mod_common
      use mod_structure

      implicit none

      INTEGER :: l, leaf, m

      REAL(kind=8) :: dvol
      REAL(kind=8), DIMENSION(nvar) :: du

      REAL(kind=8) :: pression, p_ndt, rho_ec, fac
      REAL(kind=8) :: vorticite, gpgrho
      REAL(kind=8), DIMENSION(ndim) :: grad_rho, grad_P, grad_T
      REAL(kind=8), DIMENSION(ndim,ndim) :: grad_u

      TYPE(structure_maille), pointer :: maille

! *****On parcourt les feuilles reelles (non-fictives)

      Do leaf = 1, compteur_feuille

        maille => feuille_reelle(leaf)%ptr

! +++++Controle de la densite > 0

        if( maille%u(1) < 0. ) then

          print*, ' rho < 0 : ',maille%maille_niv, &
                 maille%maille_i

          du(1) = maille%u(1) - maille%u_ndt(1)
          maille%u(1) = max ( maille%u(1), 0. )

          fac = (maille%u(1)-maille%u_ndt(1)) / du(1)

          do m = 2, nvar
            maille%u(m) = maille%u_ndt(m) + &
                          (maille%u(m) - maille%u_ndt(m))*fac
          enddo

          STOP
        endif

! +++++Controle de la pression > 0

        rho_ec = 0.
        do m = 2, nvar-1
          rho_ec = rho_ec + maille%u(m)**2
        enddo
        rho_ec = 0.5 * rho_ec / maille%u(1)

        pression = gamma*(gamma-1.) * ( maille%u(nvar) - rho_ec ) * &
                   Mach**2

        if( pression < 0. ) then

          print*, ' Pression < 0 : ',maille%maille_niv, &
                 maille%maille_i

          rho_ec = 0.
          do m = 2, nvar-1
            rho_ec = rho_ec + maille%u_ndt(m)**2
          enddo
          rho_ec = 0.5 * rho_ec / maille%u_ndt(1)

          p_ndt = gamma*(gamma-1.) * ( maille%u_ndt(nvar) - rho_ec ) * &
                  Mach**2

          fac = -.5 * p_ndt / (pression-p_ndt)

          do m = 1, nvar
            maille%u(m) = maille%u_ndt(m) + &
                          (maille%u(m) - maille%u_ndt(m))*fac
          enddo

!!          STOP
        endif

! +++++calcul des residus

        do m = 1, nvar
          du(m) = maille%u(m) - maille%u_ndt(m)
        enddo

! +++++Volume de la maille

        dvol = 1.
        do l = 1, ndim
          dvol = dvol * maille%dx(l)
        enddo

! +++++Volume total
        Volume = Volume + dvol

! +++++calcul norme L2 des residus

        do m = 1, nvar

          if( abs(du(m)) >= err_max(m) ) then
            err_max(m) = abs(du(m))
            ijk_resl2max(:, m) = maille%maille_i(:)
          endif

!!          errol2(m) = errol2(m) + (du(m)/dt)**2
          errol2(m) = errol2(m) + du(m)**2

        enddo

! *****evaluation des integrales pour interaction spot-choc

        if( ndim == 2 ) then
!           call gradient(maille, grad_rho, grad_u, grad_P, grad_T)
!
!           vorticite = grad_u(2,1) - grad_u(1,2)
!           gpgrho = ( grad_rho(1)*grad_P(2) - grad_rho(2)*grad_P(1) ) &
!                    / maille%u(1)**2
!
!           integral(1) = integral(1) +  vorticite * dvol
!
!           integral(2) = integral(2) +  abs(vorticite) * dvol
!
!           integral(3) = integral(3) + gpgrho/(gamma*Mach**2) * dvol
!
!           integral(4) = integral(4) + abs(gpgrho)/(gamma*Mach**2) * dvol
          integral(:) = integral(:) + maille%u(:) * dvol
!!
!!!!          integral(3) = max( integral(3) , maille%u(3)/maille%u(1) )
!!!!          integral(4) = max( integral(4) , pression/(gamma*Mach**2) )
!!
        else

          integral(:) = integral(:) + maille%u(:) * dvol

        endif

      Enddo

      end subroutine residu
