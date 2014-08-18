      subroutine integration_euler_acc(dir,dtt)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: dir

      INTEGER :: l, l_dt, leaf, m
      INTEGER, DIMENSION(4,3) :: ijk_dt

      REAL(kind=8) :: sigma, rho_ec, dtt
      REAL(kind=8), dimension(nvar) :: wn,wnp1

      REAL(kind=8), DIMENSION(2) :: fmp_gauche, fmp_droit,f_gauche,f_droit

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

! +++++Transfert de u(n) ===> w(n)

			wn(1) = maille%u(1)
			do m = 2,nvar
				wn(m) = maille%u(m)/wn(1)
			end do

! =====Flux Euler + correction MP a droite

			do l = 1,2
				fmp_droit(l) = maille%fmp_dt(l)
				f_droit(l) = maille%f_dt(l)

! =====Flux Euler + correction MP a gauche

! .....par defaut, on garde le flux a gauche
				fmp_gauche(l) = maille%fmp_gch(l)
				f_gauche(l) = maille%f_gch(l)
			end do

! .....Si la feuille precedente existe
        if( associated(maille%support(dir)%precedent) ) then

          local => maille%support(dir)%precedent

! .....si la feuille precedente est reelle
          if( local%feuille ) then

! .....on transfert le flux a droite pour conservation
						do l = 1,2
							fmp_gauche(l) = local%fmp_dt(l)
							f_gauche(l) = local%f_dt(l)
						end do
          endif

! .....Si la feuille precedente est a un niveau plus eleve

          if( associated(local%fils(1)%ptr) ) then

! .....si les fils de la feuille precedente sont reels
            if( local%fils(1)%ptr%feuille ) then

! .....Somme des flux Euler pour conservation

				do l = 1,2
              fmp_gauche(l) = 0.
              f_gauche(l) = 0.
            end do

              do l = 1, 2**(ndim-1)

					do m = 1,2
                fmp_gauche(m) = fmp_gauche(m) + &
                                local%fils( ijk_dt(l,dir) )%ptr%fmp_dt(m)
                f_gauche(m) = f_gauche(m) + &
                                local%fils( ijk_dt(l,dir) )%ptr%f_dt(m)
					end do

              enddo

				do l = 1,2
              fmp_gauche(l) = fmp_gauche(l) / float(2**(ndim-1))
              f_gauche(l) = f_gauche(l) / float(2**(ndim-1))
              end do

            endif

          endif

        else

          if( maille%maille_i(dir) /= 1 ) then
            Print*, ' Integration - Pb de CAL < : maille niv= ', &
             maille%maille_niv, ' i = ', maille%maille_i
            Stop

          else

! *****Conditions a la limite amont dans la direction dir

!             if( caltype_deb(dir) /= 2 .and. caltype_deb(dir) /= -1 ) &
!               call caldeb_flux_euler_ac(maille, dir, fmp_gauche,f_gauche)

! !               print*, maille%maille_i(:), f_gauche(1:2), f_droit(1:2), fmp_gauche(1:2), fmp_droit(1:2)

				endif

        endif

! *****Conditions a la limite aval dans la direction dir

        if( maille%maille_i(dir) == &
          2**maille%maille_niv*nbre_racine(dir) ) then

!           if( caltype_fin(dir) /= 2 .and. caltype_fin(dir) /= -1 ) &
!             call calfin_flux_euler_ac(maille, dir, fmp_droit,f_droit)

! !             print*,  maille%maille_i(:),f_gauche(1:2), f_droit(1:2), fmp_gauche(1:2), fmp_droit(1:2)

        endif

! *****Correction OS7 au schema Lax-Wendroff

        wnp1(1) = wn(1) - sigma*( f_droit(1) + f_gauche(1) + fmp_droit(1) - fmp_gauche(1) )
        do m = 2,nvar
			wnp1(m) = wn(m)
        end do
        wnp1(dir+1) = wn(dir+1) - sigma*( f_droit(2) + f_gauche(2) + fmp_droit(2) - fmp_gauche(2) )

! +++++Transfert de wnp1 ===> unp1

			maille%unp1(1) = wnp1(1)
			do m = 2,nvar
				maille%unp1(m) = wnp1(m)*wnp1(1)
			end do

! +++++Traitement special pour la face d'entree

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

      end subroutine integration_euler_acc
