      subroutine integration_euler(dir,dtt)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: dir

      INTEGER :: l, l_dt, leaf, m
      INTEGER, DIMENSION(4,3) :: ijk_dt

      REAL(kind=8) :: sigma, rho_ec, dtt

      REAL(kind=8), DIMENSION(nvar) :: fmp_gauche, fmp_droit,f_gauche,f_droit

      TYPE(structure_maille), pointer :: maille
      TYPE(structure_maille), pointer ::  local

! .....definition des indices des cellules a droite :
!      ijk_dt( 2**(ndim-1) , dir )

      ijk_dt(1,1) = 2
      ijk_dt(2,1) = 4
      ijk_dt(3,1) = 6
      ijk_dt(4,1) = 8

      ijk_dt(1,2) = 3
      ijk_dt(2,2) = 4
      ijk_dt(3,2) = 7
      ijk_dt(4,2) = 8

      ijk_dt(1,3) = 5
      ijk_dt(2,3) = 6
      ijk_dt(3,3) = 7
      ijk_dt(4,3) = 8

! *****On parcourt les feuilles reelles (non-fictives)

      Do leaf = 1, compteur_feuille

        maille => feuille_reelle(leaf)%ptr

! ============ Integration Euler + Correction Lax-Wendroff  ========

! .....Rapport : sigma = dt_ac / dx
        sigma = dtt / maille%dx(dir)

! +++++Maille courante

        local => maille

! =====Flux Euler + correction MP a droite

        fmp_droit(:) = maille%fmp_dt(:)
        f_droit(:) = maille%f_dt(:)

! =====Flux Euler + correction MP a gauche

! .....par defaut, on garde le flux a gauche
        fmp_gauche(:) = maille%fmp_gch(:)
        f_gauche(:) = maille%f_gch(:)

! .....Si la feuille precedente existe
        if( associated(maille%support(dir)%precedent) ) then

          local => maille%support(dir)%precedent

! .....si la feuille precedente est reelle
          if( local%feuille ) then

! .....on transfert le flux a droite pour conservation
            fmp_gauche(:) = local%fmp_dt(:)
            f_gauche(:) = local%f_dt(:)

          endif

! .....Si la feuille precedente est a un niveau plus eleve

          if( associated(local%fils(1)%ptr) ) then

! .....si les fils de la feuille precedente sont reels
            if( local%fils(1)%ptr%feuille ) then

! .....Somme des flux Euler pour conservation

              fmp_gauche(:) = 0.
              f_gauche(:) = 0.

              do l = 1, 2**(ndim-1)

                fmp_gauche(:) = fmp_gauche(:) + &
                                local%fils( ijk_dt(l,dir) )%ptr%fmp_dt(:)
                f_gauche(:) = f_gauche(:) + &
                                local%fils( ijk_dt(l,dir) )%ptr%f_dt(:)

              enddo

              fmp_gauche(:) = fmp_gauche(:) / float(2**(ndim-1))
              f_gauche(:) = f_gauche(:) / float(2**(ndim-1))

            endif

          endif

        else

          if( maille%maille_i(dir) /= 1 ) then
            Print*, ' Integration - Pb de CAL < : maille niv= ', &
             maille%maille_niv, ' i = ', maille%maille_i
            Stop

          else

! *****Conditions a la limite amont dans la direction dir

            if( caltype_deb(dir) /= 2 .and. caltype_deb(dir) /= -1 ) &
              call caldeb_flux_euler(maille, dir, fmp_gauche)

          endif

        endif

! *****Conditions a la limite aval dans la direction dir

        if( maille%maille_i(dir) == &
          2**maille%maille_niv*nbre_racine(dir) ) then

          if( caltype_fin(dir) /= 2 .and. caltype_fin(dir) /= -1 ) &
            call calfin_flux_euler(maille, dir, fmp_droit)

        endif

! *****Correction OSMP7 au schema Lax-Wendroff

        maille%unp1(:) = maille%u(:) - sigma*( f_droit(:) + f_gauche(:) + fmp_droit(:) - fmp_gauche(:) )

! +++++Traitement special pour la cavité entrainée

#if Cav1 .or. Cav2
        if( caltype_fin(2) == 0 .and. maille%maille_i(2) == 2**maille%maille_niv*nbre_racine(2) ) then

          maille%unp1(1) = rho_right
          rho_ec = 0.
          do m = 2,nvar-1
            maille%unp1(m) = rho_right*u_cav(maille%maille_i(1),m-1 )
            rho_ec = rho_ec + 0.5*rho_right*maille%unp1(m)**2
          enddo
          maille%unp1(nvar) = p_left/(gamma*(gamma-1.)*Mach**2) + rho_ec

        endif
#endif
      Enddo

      end subroutine integration_euler
