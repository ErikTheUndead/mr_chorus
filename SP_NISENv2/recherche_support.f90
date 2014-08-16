      subroutine recherche_support(maille, s_st, f_support)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      TYPE(structure_maille), POINTER :: maille

      INTEGER :: i, j, k, l, lm, m, n
      INTEGER :: i_sym, j_sym, k_sym
      INTEGER :: indice, increment
      INTEGER :: niv_sup, s_st, sj, sk
      INTEGER :: si_deb, sj_deb, sk_deb, si_fin, sj_fin, sk_fin
      INTEGER, DIMENSION(ndim) :: ind_loc, ind_sup, ind_cor
      INTEGER, DIMENSION(3) :: ind

      REAL(kind=8) :: rho_ec, rho_wall, rhoE_wall

      REAL(kind=8), DIMENSION(-s_st:s_st,-s_st:s_st,-s_st:s_st,nvar) :: &
                    f_support

      TYPE(structure_maille), pointer :: courant

! =======================================================================

! .....niveau et indice de la maille courante

      niv_sup = maille%maille_niv
      ind_sup(:) = maille%maille_i(:)

! +++++Definition des bornes debut et fin

      if( caltype_deb(1) == 2 ) then
        si_deb = -s_st
        si_fin =  s_st
      else 
        si_deb = max(1, ind_sup(1)-s_st) - ind_sup(1)
        si_fin = min(2**niv_sup*nbre_racine(1), ind_sup(1)+s_st) - ind_sup(1)
      endif

      if(ndim == 3) then
        sj = s_st
        sk = s_st

        if( caltype_deb(2) == 2 ) then
          sj_deb = -sj
          sj_fin =  sj
        else
          sj_deb = max(1, ind_sup(2)-sj) - ind_sup(2)
          sj_fin = min(2**niv_sup*nbre_racine(2), ind_sup(2)+sj) - ind_sup(2)
        endif

        if( caltype_deb(3) == 2 ) then
          sk_deb = -sk
          sk_fin =  sk
        else
          sk_deb = max(1, ind_sup(3)-sk) - ind_sup(3)
          sk_fin = min(2**niv_sup*nbre_racine(3), ind_sup(3)+sk) - ind_sup(3)
        endif

      elseif(ndim == 2) then
        sj = s_st
        sk = 0

        if( caltype_deb(2) == 2 ) then
          sj_deb = -sj
          sj_fin =  sj
        else
          sj_deb = max(1, ind_sup(2)-sj) - ind_sup(2)
          sj_fin = min(2**niv_sup*nbre_racine(2), ind_sup(2)+sj) - ind_sup(2)
        endif

        sk_deb = 0
        sk_fin = 0

      elseif(ndim == 1) then
        sj = 0
        sk = 0

        sj_deb = 0
        sj_fin = 0

        sk_deb = 0
        sk_fin = 0

      endif

! ======== Recherche d'un support d'interpolation [-s_st , +s_st] ========

! +++++initialisation
!!      f_support(:,:,:,:) = 0.

! .....Valorisation pour la maille courante
      do m = 1,nvar
        f_support(0,0,0,m) = maille%u(m)
      enddo
          
! *****Recherche de tous les points du support et valorisation 

! .....polynome de degre >= 1
      if(s_st > 0 ) then

! +++++Recherche des points dan sle domaine de calcul

        do k = sk_deb, sk_fin
          ind(3) =  k

          do j = sj_deb, sj_fin
            ind(2) = j

            do i = si_deb, si_fin
              ind(1) = i

              do l = 1,ndim
                ind_loc(l) = ind_sup(l) + ind(l)

! .....Cas de la periodicite
                if( caltype_deb(l) == 2 ) then

                  do while(ind_loc(l) <= 0 .or. &
                    ind_loc(l) > 2**niv_sup*nbre_racine(l) )
!
                    if(ind_loc(l) <= 0 ) &
                      ind_cor(l) = 2**niv_sup*nbre_racine(l) + ind_loc(l)
                    if(ind_loc(l) > 2**niv_sup ) &
                      ind_cor(l) = ind_loc(l) - 2**niv_sup*nbre_racine(l)
!
                    ind_loc(l) = ind_cor(l)
                  enddo
                endif

              enddo

! .....recherche du point dans l'arbre 
              indice = 1
              increment = 1
              do l = 1, ndim
                indice = indice + (ind_loc(l) - 1) * increment
                increment = increment * 2**niv_sup * nbre_racine(l)
              enddo
              indice = indice + (2**(ndim*niv_sup) - 1)/(2**ndim-1) * &
                                nbre_racine(1)*nbre_racine(2)*nbre_racine(3)

              if( .not.associated(hash_table(indice)%ptr) ) then

                print*, ' !! Recherche Support !! : ', &
                        ' Pb de recherche support maille - ', &
                        'Maille niv = ',niv_sup, ' Reelle = ',maille%feuille, &
                        ' Maille ijk = ',ind_sup, &
                        ' i, j, k loc = ',ind_loc

                stop
              endif

              courant => hash_table(indice)%ptr          

              do m = 1,nvar
                f_support(i, j, k, m) = courant%u(m)
              enddo

            enddo
          enddo
        enddo

! ========================= CAL suivant X ==========================

! +++++Conditions aux limites amont

        if( si_deb > -s_st ) then

! .....Conditions de symetrie

          if( caltype_deb(1) == -1 ) then

            do k = -sk, sk
              do j = -sj, sj

                do i = -s_st, si_deb-1
                  i_sym = 2*si_deb-(i+1)
                  do m = 1, nvar
                    f_support(i,j,k,m) = f_support(i_sym,j,k,m)
                  enddo
                  f_support(i,j,k,2) = -f_support(i_sym,j,k,2)
                enddo

              enddo
            enddo

! .....Conditions de frontiere fluide ==> derivee nulle

          elseif( caltype_deb(1) == 0 ) then

            do k = -sk, sk
              do j = -sj, sj

                do i = -s_st, si_deb-1
                  do m = 1, nvar
                    f_support(i,j,k,m) = f_support(si_deb,j,k,m)
                  enddo
                enddo

              enddo
            enddo

! .....Conditions de paroi

          elseif( caltype_deb(1) == 1 ) then

            if( cle_twall ) then
              do k = -sk, sk
                do j = -sj, sj

! .....paroi adiabatique ==> gradient densite et presssion  = 0
                  rho_wall = f_support(si_deb,j,k,1)
                  rho_ec = 0.
                  do m = 2, nvar-1
                    rho_ec = rho_ec + f_support(si_deb,j,k,m)**2               
                  enddo
                  rhoE_wall = f_support(si_deb,j,k,nvar) - &
                              0.5*rho_ec/f_support(si_deb,j,k,1)

                  if( si_deb+1 <= si_fin ) then
                    rho_wall = (9.*rho_wall - &
                                   f_support(si_deb+1,j,k,1))/8.
                    rho_ec = 0.
                    do m = 2, nvar-1
                      rho_ec = rho_ec + f_support(si_deb+1,j,k,m)**2
                    enddo
                    rhoE_wall = (9.*rhoE_wall  - &
                               (f_support(si_deb+1,j,k,nvar) - &
                                0.5*rho_ec/f_support(si_deb+1,j,k,1)) )/8.
                  endif

                  do i = -s_st, si_deb-1
                    i_sym = 2*si_deb-(i+1)
                    f_support(i,j,k,1) = rho_wall
                    do m = 2, nvar-1
!!                      f_support(i,j,k,m) = 0.
                      f_support(i,j,k,m) = - f_support(i_sym,j,k,m)
                    enddo
                    f_support(i,j,k,nvar) = rhoE_wall
                  enddo
            
                enddo
              enddo
            else
! .....paroi  thermostatee 
              do k = -sk, sk
                do j = -sj, sj

                  do i = -s_st, si_deb-1
                    i_sym = 2*si_deb-(i+1)
                    f_support(i,j,k,1) = f_support(i_sym,j,k,1)
                    do m = 2, nvar-1
                      f_support(i,j,k,m) = - f_support(i_sym,j,k,m)
                    enddo
                    f_support(i,j,k,nvar) = f_support(i_sym,j,k,nvar)
                  enddo

                enddo
              enddo

            endif
          endif

        endif

! *****Conditions a la limite aval

        if( si_fin < s_st ) then

! .....Conditions de symetrie

          if( caltype_fin(1) == -1 ) then

            do k = -sk, sk
              do j = -sj, sj

                do i = si_fin+1, s_st
                  i_sym = 2*si_fin-(i-1)
                  do m = 1, nvar
                    f_support(i,j,k,m) = f_support(i_sym,j,k,m)
                  enddo
                  f_support(i,j,k,2) = -f_support(i_sym,j,k,2)
                enddo

              enddo
            enddo

! .....Conditions de frontiere fluide ==> derivee nulle
          elseif( caltype_fin(1) == 0 ) then

            do k = -sk, sk
              do j = -sj, sj

                do i = si_fin+1, s_st
                  do m = 1, nvar
                    f_support(i,j,k,m) = f_support(si_fin,j,k,m)
                  enddo
                enddo

              enddo
            enddo

! .....Conditions de paroi
          elseif( caltype_fin(1) == 1 ) then

            if( cle_twall ) then
              do k = -sk, sk
                do j = -sj, sj

! .....paroi adiabatique ==> gradient densite et presssion  = 0
                  rho_wall = (9.*f_support(si_fin,j,k,1) - &
                                 f_support(si_fin-1,j,k,1))/8.
                  rho_ec = 0.
                  do m = 2, nvar-1
                    rho_ec = rho_ec + f_support(si_fin,j,k,m)**2               
                  enddo
                  rhoE_wall = 9.*(f_support(si_fin,j,k,nvar) - &
                                  0.5*rho_ec/f_support(si_fin,j,k,1))
                  rho_ec = 0.
                  do m = 2, nvar-1
                    rho_ec = rho_ec + f_support(si_fin-1,j,k,m)**2
                  enddo
                  rhoE_wall = (rhoE_wall  - &
                           (f_support(si_fin-1,j,k,nvar) - &
                            0.5*rho_ec/f_support(si_fin-1,j,k,1)) )/8.
                  do i = si_fin+1, s_st
                    i_sym = 2*si_fin-(i-1)
                    f_support(i,j,k,1) = rho_wall
                    do m = 2, nvar-1
!!                      f_support(i,j,k,m) = 0.
                      f_support(i,j,k,m) = - f_support(i_sym,j,k,m)
                    enddo
                    f_support(i,j,k,nvar) = rhoE_wall
                  enddo
                enddo
              enddo
            
            else
! .....paroi  thermostatee 
              do k= -sk, sk
                do j = -sj, sj

                  do i = si_fin+1, s_st
                    i_sym = 2*si_fin-(i-1)
                    f_support(i,j,k,1) = f_support(i_sym,j,k,1)
                    do m = 2, nvar-1
                      f_support(i,j,k,m) = - f_support(i_sym,j,k,m)
                    enddo
                    f_support(i,j,k,nvar) = f_support(i_sym,j,k,nvar)
                  enddo
                enddo
              enddo

            endif

          endif

        endif

! ========================= CAL suivant Y ==========================

! +++++Conditions aux limites amont

        if( ndim >= 2 .and. sj_deb > -sj ) then

! .....Conditions de symetrie

          if( caltype_deb(2) == -1 ) then

            do k = -sk, sk
              do i = -s_st, s_st

                do j = -sj, sj_deb-1
                  j_sym = 2*sj_deb-(j+1)
                  do m = 1, nvar
                    f_support(i,j,k,m) = f_support(i,j_sym,k,m)
                  enddo
                  f_support(i,j,k,3) = -f_support(i,j_sym,k,3)
                enddo

              enddo
            enddo

! .....Conditions de frontiere fluide ==> derivee nulle

          elseif( caltype_deb(2) == 0 ) then

            do k = -sk, sk
              do i = -s_st, s_st

                do j = -sj, sj_deb-1
                  do m = 1, nvar
                    f_support(i,j,k,m) = f_support(i,sj_deb,k,m)
                  enddo
                enddo

              enddo
            enddo

! .....Conditions de paroi

          elseif( caltype_deb(2) == 1 ) then

            if( cle_twall ) then
              do k = -sk, sk
                do i = -s_st, s_st

! .....paroi adiabatique ==> gradient densite et presssion  = 0
                  rho_wall = f_support(i,sj_deb,k,1)

                  rho_ec = 0.
                  do m = 2, nvar-1
                    rho_ec = rho_ec + f_support(i,sj_deb,k,m)**2               
                  enddo
                  rhoE_wall = f_support(i,sj_deb,k,nvar) - &
                              0.5*rho_ec/f_support(i,sj_deb,k,1)

                  if( sj_deb+1 <= sj_fin ) then
                    rho_wall = (9.*rho_wall - f_support(i,sj_deb+1,k,1))/8.
                    rho_ec = 0.
                    do m = 2, nvar-1
                      rho_ec = rho_ec + f_support(i,sj_deb+1,k,m)**2 
                    enddo
                    rhoE_wall = (9.*rhoE_wall  - &
                                 (f_support(i,sj_deb+1,k,nvar) - &
                                  0.5*rho_ec/f_support(i,sj_deb+1,k,1)) )/8.
                  endif

                  do j = -sj, sj_deb-1
                    j_sym = 2*sj_deb-(j+1)
                    f_support(i,j,k,1) = rho_wall
                    do m = 2, nvar-1
!!                      f_support(i,j,k,m) = 0.
                      f_support(i,j,k,m) = - f_support(i,j_sym,k,m)
                    enddo
                    f_support(i,j,k,nvar) = rhoE_wall
                  enddo
            
                enddo
              enddo
            else
! .....paroi  thermostatee 
              do k = -sk, sk
                do i = -s_st, s_st

                  do j = -sj, sj_deb-1
                    j_sym = 2*sj_deb-(j+1)
                    f_support(i,j,k,1) = f_support(i,j_sym,k,1)
                    do m = 2, nvar-1
                      f_support(i,j,k,m) = - f_support(i,j_sym,k,m)
                    enddo
                    f_support(i,j,k,nvar) = f_support(i,j_sym,k,nvar)
                  enddo

                enddo
              enddo

            endif
          endif

        endif

! *****Conditions a la limite aval

        if( ndim >= 2 .and. sj_fin < sj ) then

! .....Conditions de symetrie

          if( caltype_fin(2) == -1 ) then

            do k = -sk, sk
              do i = -s_st, s_st

                do j = sj_fin+1, sj
                  j_sym = 2*sj_fin-(j-1)
                  do m = 1, nvar
                    f_support(i,j,k,m) = f_support(i,j_sym,k,m)
                  enddo
                  f_support(i,j,k,3) = -f_support(i,j_sym,k,3)
                enddo

              enddo
            enddo

! .....Conditions de frontiere fluide ==> derivee nulle
          elseif( caltype_fin(2) == 0 ) then

            do k = -sk, sk
              do i = -s_st, s_st

                do j = sj_fin+1, sj
                  do m = 1, nvar
                    f_support(i,j,k,m) = f_support(i,sj_fin,k,m)
                  enddo
                enddo

              enddo
            enddo

! .....Conditions de paroi
          elseif( caltype_fin(2) == 1 ) then

            if( cle_twall ) then
              do k = -sk, sk
                do i = -s_st, s_st

! .....paroi adiabatique ==> gradient densite et presssion  = 0
                  rho_wall = (9.*f_support(i,sj_fin,k,1) - &
                                 f_support(i,sj_fin-1,k,1))/8.
                  rho_ec = 0.
                  do m = 2, nvar-1
                    rho_ec = rho_ec + f_support(i,sj_fin,k,m)**2               
                  enddo
                  rhoE_wall = 9.*(f_support(i,sj_fin,k,nvar) - &
                                  0.5*rho_ec/f_support(i,sj_fin,k,1))
                  rho_ec = 0.
                  do m = 2, nvar-1
                    rho_ec = rho_ec + f_support(i,sj_fin-1,k,m)**2
                  enddo
                  rhoE_wall = (rhoE_wall  - &
                           (f_support(i,sj_fin-1,k,nvar) - &
                            0.5*rho_ec/f_support(i,sj_fin-1,k,1)) )/8.
                  do j = sj_fin+1, sj
                    j_sym = 2*sj_fin-(j-1)
                    f_support(i,j,k,1) = rho_wall
                    do m = 2, nvar-1
!!                      f_support(i,j,k,m) = 0.
                      f_support(i,j,k,m) = - f_support(i,j_sym,k,m)
                    enddo
                    f_support(i,j,k,nvar) = rhoE_wall
                  enddo
                enddo
              enddo
            
            else
! .....paroi  thermostatee 
              do k= -sk, sk
                do i = -s_st, s_st

                  do j = sj_fin+1, sj
                    j_sym = 2*sj_fin-(j-1)
                    f_support(i,j,k,1) = f_support(i,j_sym,k,1)
                    do m = 2, nvar-1
                      f_support(i,j,k,m) = - f_support(i,j_sym,k,m)
                    enddo
                    f_support(i,j,k,nvar) = f_support(i,j_sym,k,nvar)
                  enddo
                enddo
              enddo

            endif

          endif

        endif

! ========================= CAL suivant Z ==========================

! +++++Conditions aux limites amont

        if( ndim >= 3 .and. sk_deb > -sk ) then

! .....Conditions de symetrie

          if( caltype_deb(3) == -1 ) then

            do j = -sj, sj
              do i = -s_st, s_st

                do k = -sk, sk_deb-1
                  k_sym = 2*sk_deb-(k+1)
                  do m = 1, nvar
                    f_support(i,j,k,m) = f_support(i,j,k_sym,m)
                  enddo
                  f_support(i,j,k,4) = -f_support(i,j,k_sym,4)
                enddo

              enddo
            enddo

! .....Conditions de frontiere fluide ==> derivee nulle

          elseif( caltype_deb(3) == 0 ) then

            do j = -sj, sj
              do i = -s_st, s_st

                do k = -sk, sk_deb-1
                  do m = 1, nvar
                    f_support(i,j,k,m) = f_support(i,j,sk_deb,m)
                  enddo
                enddo

              enddo
            enddo

! .....Conditions de paroi

          elseif( caltype_deb(3) == 1 ) then

            if( cle_twall ) then
              do j = -sj, sj
                do i = -s_st, s_st

! .....paroi adiabatique ==> gradient densite et presssion  = 0
                  rho_wall = f_support(i,j,sk_deb,1)

                  rho_ec = 0.
                  do m = 2, nvar-1
                    rho_ec = rho_ec + f_support(i,j,sk_deb,m)**2               
                  enddo
                  rhoE_wall = f_support(i,j,sk_deb,nvar) - &
                              0.5*rho_ec/f_support(i,j,sk_deb,1)

                  if( sk_deb+1 <= sk_fin ) then
                    rho_wall = (9.*rho_wall - f_support(i,j,sk_deb+1,1))/8.
                    rho_ec = 0.
                    do m = 2, nvar-1
                      rho_ec = rho_ec + f_support(i,j,sk_deb+1,m)**2
                    enddo
                    rhoE_wall = (9.*rhoE_wall  - &
                                 (f_support(i,j,sk_deb+1,nvar) - &
                                  0.5*rho_ec/f_support(i,j,sk_deb+1,1)) )/8.
                  endif

                  do k = -sk, sk_deb-1
                    k_sym = 2*sk_deb-(k+1)
                    f_support(i,j,k,1) = rho_wall
                    do m = 2, nvar-1
!!                      f_support(i,j,k,m) = 0.
                      f_support(i,j,k,m) = - f_support(i,j,k_sym,m)
                    enddo
                    f_support(i,j,k,nvar) = rhoE_wall
                  enddo
            
                enddo
              enddo
            else
! .....paroi  thermostatee 
              do j = -sj, sj
                do i = -s_st, s_st

                  do k = -sk, sk_deb-1
                    k_sym = 2*sk_deb-(k+1)
                    f_support(i,j,k,1) = f_support(i,j,k_sym,1)
                    do m = 2, nvar-1
                      f_support(i,j,k,m) = - f_support(i,j,k_sym,m)
                    enddo
                    f_support(i,j,k,nvar) = f_support(i,j,k_sym,nvar)
                  enddo

                enddo
              enddo

            endif
          endif

        endif

! *****Conditions a la limite aval

        if( ndim >= 3 .and. sk_fin < sk ) then

! .....Conditions de symetrie

          if( caltype_fin(3) == -1 ) then

            do j = -sj, sj
              do i = -s_st, s_st

                do k = sk_fin+1, sk
                  k_sym = 2*sk_fin-(k-1)
                  do m = 1, nvar
                    f_support(i,j,k,m) = f_support(i,j,k_sym,m)
                  enddo
                  f_support(i,j,k,4) = -f_support(i,j,k_sym,4)
                enddo

              enddo
            enddo

! .....Conditions de frontiere fluide ==> derivee nulle
          elseif( caltype_fin(3) == 0 ) then

            do j = -sj, sj
              do i = -s_st, s_st

                do k = sk_fin+1, sk
                  do m = 1, nvar
                    f_support(i,j,k,m) = f_support(i,j,sk_fin,m)
                  enddo
                enddo

              enddo
            enddo

! .....Conditions de paroi
          elseif( caltype_fin(3) == 1 ) then

            if( cle_twall ) then
              do j = -sj, sj
                do i = -s_st, s_st

! .....paroi adiabatique ==> gradient densite et presssion  = 0
                  rho_wall = (9.*f_support(i,j,sk_fin,1) - &
                                 f_support(i,j,sk_fin-1,1))/8.
                  rho_ec = 0.
                  do m = 2, nvar-1
                    rho_ec = rho_ec + f_support(i,j,sk_fin,m)**2
                  enddo
                  rhoE_wall = 9.*(f_support(i,j,sk_fin,nvar) - &
                                  0.5*rho_ec/f_support(i,j,sk_fin,1))
                  rho_ec = 0.
                  do m = 2, nvar-1
                    rho_ec = rho_ec + f_support(i,j,sk_fin-1,m)**2
                  enddo
                  rhoE_wall = (rhoE_wall  - &
                           (f_support(i,j,sk_fin-1,nvar) - &
                            0.5*rho_ec/f_support(i,j,sk_fin-1,1)) )/8.
                  do k = sk_fin+1, sk
                    k_sym = 2*sk_fin-(k-1)
                    f_support(i,j,k,1) = rho_wall
                    do m = 2, nvar-1
!!                      f_support(i,j,k,m) = 0.
                      f_support(i,j,k,m) = - f_support(i,j,k_sym,m)
                    enddo
                    f_support(i,j,k,nvar) = rhoE_wall
                  enddo
                enddo
              enddo
            
            else
! .....paroi  thermostatee 
              do j= -sj, sj
                do i = -s_st, s_st

                  do k = sk_fin+1, sk
                    k_sym = 2*sk_fin-(k-1)
                    f_support(i,j,k,1) = f_support(i,j,k_sym,1)
                    do m = 2, nvar-1
                      f_support(i,j,k,m) = - f_support(i,j,k_sym,m)
                    enddo
                    f_support(i,j,k,nvar) = f_support(i,j,k_sym,nvar)
                  enddo
                enddo
              enddo

            endif

          endif

        endif

      endif

      end subroutine recherche_support
