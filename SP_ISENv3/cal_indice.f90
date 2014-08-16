	subroutine cal_indice(niv_loc, ind_loc, ind)
	
        use mod_common

	IMPLICIT NONE
	
        INTEGER :: l
        INTEGER :: niv_loc
        INTEGER, DIMENSION(ndim) :: ind, ind_loc

! =========== gestion des indices pour une recherche dans l'arbre =========

! +++++Application des conditions aux limites pour
!      un indice (ind_loc) a un niveau (niv_loc) donnes

        ind(:) = ind_loc(:)

        do l = 1, ndim

! .....conditions aux limites ==> symetrie
          if( caltype_deb(l) <= 1 .OR. caltype_fin(l) <= 1 ) then

            do while(ind_loc(l) <= 0 .or. &
              ind_loc(l) > 2**niv_loc*nbre_racine(l) )

              if(ind_loc(l) <= 0 ) ind(l) = 1 - ind_loc(l)
              if(ind_loc(l) > 2**niv_loc*nbre_racine(l) ) &
                ind(l) = 2**(niv_loc+1)*nbre_racine(l) - ind_loc(l) + 1

              ind_loc(l) = ind(l)
            enddo

! .....conditions aux limites ==> periodicite
          elseif( caltype_deb(l) == 2 .OR. caltype_fin(l) == 2 ) then

            do while(ind_loc(l) <= 0 .or. &
              ind_loc(l) > 2**niv_loc*nbre_racine(l) )

              if(ind_loc(l) <= 0 ) &
                ind(l) = 2**niv_loc*nbre_racine(l) + ind_loc(l)
              if(ind_loc(l) > 2**niv_loc ) &
                ind(l) = ind_loc(l) - 2**niv_loc*nbre_racine(l)

              ind_loc(l) = ind(l)
            enddo

          endif
        enddo

	end subroutine cal_indice
