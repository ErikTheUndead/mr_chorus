      subroutine flux_visc(dir)

      use mod_common
      use mod_structure
      use mod_fonction

      IMPLICIT NONE

      INTEGER :: dir

      LOGICAL :: flag
 
      INTEGER :: i_dir, i_sym, l, l_dt, leaf, m
      INTEGER :: s_deb, s_fin
      INTEGER :: ijk_gch
      INTEGER, DIMENSION(3) :: i_beg, i_end, i_pas
      INTEGER, DIMENSION(4,3) :: ijk_dt

      REAL(kind=8), DIMENSION(nvar, -1:1, ndim) :: u_0
      REAL(kind=8), DIMENSION(nvar, -1:1, ndim) :: u_1

      REAL(kind=8) :: mu_0, mu_1
      REAL(kind=8) :: ec, T_0, T_1

      REAL(kind=8) :: mu_ip12
      REAL(kind=8) :: divV, T_xyz
      REAL(kind=8), DIMENSION(ndim, ndim) :: u_xyz
      REAL(kind=8), DIMENSION(ndim) :: u_ip12

      REAL(kind=8), DIMENSION(nvar) :: fvisc_droit, C_visc

      TYPE(structure_maille), pointer :: maille

      TYPE(structure_maille), pointer :: courant, crt_suivant, local

! .....Definition des indices des cellules a droite : 
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

! .....Indices de debut et fin dans les directions transverses

      i_beg(1) = 2
      i_end(1) = ndim
      i_pas(1) = 1

      i_beg(2) = 1
      i_end(2) = ndim
      i_pas(2) = 2

      i_beg(3) = 1
      i_end(3) = 2
      i_pas(3) = 1

! *****On parcourt les feuilles reelles (non-fictives)

      Do leaf = 1, compteur_feuille

        maille => feuille_reelle(leaf)%ptr

! .....Test pour reconstruire au niveau de la feuille ou au niveau superieur ?

        flag = .false.

        if( associated(maille%support(dir)%suivant) ) then

          if( associated(maille%support(dir)%suivant%fils(1)%ptr) ) then

            if(maille%support(dir)%suivant%fils(1)%ptr%feuille) &
              flag = .true.

          endif
        endif

! *****Si la feuille a des fils et que 
!      la feuille suivante a des fils non-fictifs

        if( flag .and. associated(maille%fils(1)%ptr) ) then

! ======== Flux visqueux a droite pour toutes les feuilles =====

          fvisc_droit(:) = 0.

! +++++La correction est evaluee au niveau superieur (Fils de la feuille)

          do l_dt = 1, 2**(ndim-1)

! *****Maille courante

            courant => maille%fils( ijk_dt(l_dt,dir) )%ptr

! .....vicosite dynamique et conductivite thermique
            call viscosite(courant, mu_0)

! .....vitesses au centre de la maille
            do l = 1, ndim
              do m = 1, nvar
                u_0(m, 0, l) = courant%u(m)
              enddo
            enddo

! *****Directions perpendiculaires

            if( ndim > 1 ) then

              do i_dir = i_beg(dir), i_end(dir), i_pas(dir)

! +++++Point precedent dans les directions perpendiculaires

                if( associated(courant%support(i_dir)%precedent) ) then

                  local => courant%support(i_dir)%precedent

                  do m = 1, nvar
                    u_0(m, -1, i_dir) = local%u(m)
                  enddo

! +++++Conditions aux limites (sauf periodiques)

                elseif( courant%maille_i(i_dir) == 1 ) then

! .....Conditions de symetrie

                  if( caltype_deb(i_dir) == -1 ) then

                    do m = 1, nvar
                      u_0(m, -1, i_dir) = u_0(m, 0, i_dir)
                    enddo
                    u_0(i_dir+1, -1, i_dir) = - u_0(i_dir+1, 0, i_dir)

! .....Conditions de frontiere fluide ==> derivee nulle
  
                  elseif( caltype_deb(i_dir) == 0 ) then

                    do m = 1, nvar
                      u_0(m, -1, i_dir) = u_0(m, 0, i_dir)
                    enddo

! .....Conditions de paroi

                  elseif( caltype_deb(i_dir) == 1 ) then

                    u_0(1, -1, i_dir) = u_0(1, 0, i_dir)
                    do m = 2, nvar-1
!!                      u_0(m, -1, i_dir) = 0
                      u_0(m, -1, i_dir) = - u_0(m, 0, i_dir)
                    enddo
                    u_0(nvar, -1, i_dir) = u_0(nvar, 0, i_dir)

                  endif

! +++++Probleme d'absence de maille voisine !

                else

                  Print*
                  Print*, ' Flux_visc : maille precedente absente !!', &
                          ' direction = ', dir, 'transverse = ',i_dir, &
                          ' niv = ', courant%maille_niv, &
                          ' indices = ', courant%maille_i, &
                          ' niv = ', local%maille_niv, &
                          ' indices = ', local%maille_i

                  Stop

                endif

! +++++Point suivant dans les directions perpendiculaires

                if( associated(courant%support(i_dir)%suivant) ) then

                  local => courant%support(i_dir)%suivant

                  do m = 1, nvar
                    u_0(m, 1, i_dir) = local%u(m)
                  enddo

! +++++Conditions aux limites

                elseif( courant%maille_i(i_dir) == &
                        2**courant%maille_niv * nbre_racine(i_dir) ) then

! .....Conditions de symetrie

                  if( caltype_fin(i_dir) == -1 ) then

                    do m = 1, nvar
                      u_0(m, 1, i_dir) = u_0(m, 0, i_dir)
                    enddo
                    u_0(i_dir+1, 1, i_dir) = - u_0(i_dir+1, 0, i_dir)

! .....Conditions de frontiere fluide ==> derivee nulle
  
                  elseif( caltype_fin(i_dir) == 0 ) then

                    do m = 1, nvar
                      u_0(m, 1, i_dir) = u_0(m, 0, i_dir)
                    enddo

! .....Conditions de paroi

                  elseif( caltype_fin(i_dir) == 1 ) then

                    u_0(1, 1, i_dir) = u_0(1, 0, i_dir)
                    do m = 2, nvar-1
!!                      u_0(m, 1, i_dir) = 0
                      u_0(m, 1, i_dir) = - u_0(m, 0, i_dir)
                    enddo
                    u_0(nvar, 1, i_dir) = u_0(nvar, 0, i_dir)

                  endif

! +++++Probleme d'absence de maille voisine !

                else

                  Print*
                  Print*, ' Flux_visc : maille suivante absente !!', &
                          ' direction = ', dir,' transverse = ',i_dir, &
                          ' niv = ', courant%maille_niv, &
                          ' indices = ', courant%maille_i, &
                          ' niv = ', local%maille_niv, &
                          ' indices = ', local%maille_i

                  Stop

                endif

              enddo

            endif

! *****Maille suivante

            if( associated(courant%support(dir)%suivant) ) then

              crt_suivant => courant%support(dir)%suivant

! .....vicosite dynamique et conductivite thermique
              call viscosite(crt_suivant, mu_1)

! .....vitesses au centre de la maille
              do l = 1, ndim
                do m = 1, nvar
                  u_1(m, 0, l) = crt_suivant%u(m)
                enddo
              enddo

! *****Directions perpendiculaires

              if( ndim > 1 ) then

                do i_dir = i_beg(dir), i_end(dir), i_pas(dir)

! +++++Point precedent dans les directions perpendiculaires

                  if( associated(crt_suivant%support(i_dir)%precedent) ) then

                    local => crt_suivant%support(i_dir)%precedent

                    do m = 1, nvar
                      u_1(m, -1, i_dir) = local%u(m)
                    enddo

! +++++Conditions aux limites

                  elseif( crt_suivant%maille_i(i_dir) == 1 ) then

! .....Conditions de symetrie

                    if( caltype_deb(i_dir) == -1 ) then

                      do m = 1, nvar
                        u_1(m, -1, i_dir) = u_1(m, 0, i_dir)
                      enddo
                      u_1(i_dir+1, -1, i_dir) = - u_1(i_dir+1, 0, i_dir)

! .....Conditions de frontiere fluide ==> derivee nulle
  
                    elseif( caltype_deb(i_dir) == 0 ) then

                      do m = 1, nvar
                        u_1(m, -1, i_dir) = u_1(m, 0, i_dir)
                      enddo

! .....Conditions de paroi

                    elseif( caltype_deb(i_dir) == 1 ) then

                      u_1(1, -1, i_dir) = u_1(1, 0, i_dir)
                      do m = 2, nvar-1
!!                        u_1(m, -1, i_dir) = 0
                        u_1(m, -1, i_dir) = - u_1(m, 0, i_dir)
                      enddo
                      u_1(nvar, -1, i_dir) = u_1(nvar, 0, i_dir)

                    endif

! +++++Probleme d'absence de maille voisine !

                  else

                    Print*
                    Print*, ' Flux_visc : maille precedente absente !!', &
                            ' direction = ', dir,' transverse = ',i_dir, &
                            ' niv = ', crt_suivant%maille_niv, &
                          ' indices = ', courant%maille_i, &
                          ' niv = ', local%maille_niv, &
                          ' indices = ', local%maille_i

                    Stop

                  endif

! +++++Point suivant dans les directions perpendiculaires

                  if( associated(crt_suivant%support(i_dir)%suivant) ) then

                    local => crt_suivant%support(i_dir)%suivant

                    do m = 1, nvar
                      u_1(m, 1, i_dir) = local%u(m)
                    enddo

! +++++Conditions aux limites

                  elseif( crt_suivant%maille_i(i_dir) == &
                          2**crt_suivant%maille_niv * nbre_racine(i_dir) ) then

! .....Conditions de symetrie

                    if( caltype_fin(i_dir) == -1 ) then

                      do m = 1, nvar
                        u_1(m, 1, i_dir) = u_1(m, 0, i_dir)
                      enddo
                      u_1(i_dir+1, 1, i_dir) = - u_1(i_dir+1, 0, i_dir)

! .....Conditions de frontiere fluide ==> derivee nulle
  
                    elseif( caltype_fin(i_dir) == 0 ) then

                      do m = 1, nvar
                        u_1(m, 1, i_dir) = u_1(m, 0, i_dir)
                      enddo

! .....Conditions de paroi

                    elseif( caltype_fin(i_dir) == 1 ) then

                      u_1(1, 1, i_dir) = u_1(1, 0, i_dir)
                      do m = 2, nvar-1
!!                        u_1(m, 1, i_dir) = 0
                        u_1(m, 1, i_dir) = - u_1(m, 0, i_dir)
                      enddo
                      u_1(nvar, 1, i_dir) = u_1(nvar, 0, i_dir)

                    endif

! +++++Probleme d'absence de maille voisine !

                  else

                    Print*
                    Print*, ' Flux_visc : maille suivante absente !!', &
                        ' direction = ', dir,' transverse = ', i_dir, &
                        ' niv = ', crt_suivant%maille_niv, &
                          ' indices = ', courant%maille_i, &
                          ' niv = ', local%maille_niv, &
                          ' indices = ', local%maille_i

                    Stop

                  endif

                enddo

              endif

! +++++Conditions aux limites suivant X

            elseif( courant%maille_i(dir) == &
                    2**courant%maille_niv * nbre_racine(dir) ) then

              mu_1 = mu_0

! .....Conditions de symetrie

              if( caltype_fin(dir) == -1 ) then

                do i_dir = 1, ndim
                  do l = -1, 1
                    do m = 1, nvar
                      u_1(m, l, i_dir) = u_0(m, l, i_dir)
                    enddo
                    u_1(dir+1, l, i_dir) = - u_0(dir+1, l, i_dir)
                  enddo
                enddo

! .....Conditions de frontiere fluide ==> derivee nulle
  
              elseif( caltype_fin(dir) == 0 ) then

                do i_dir = 1, ndim
                  do l = -1, 1
                    do m = 1, nvar
                      u_1(m, l, i_dir) = u_0(m, l, i_dir)
                    enddo
                  enddo
                enddo

! .....Conditions de paroi

              elseif( caltype_fin(dir) == 1 ) then

                do i_dir = 1, ndim
                  do l = -1, 1
                    u_1(1, l, i_dir) = u_0(1, l, i_dir)
                    do m = 2, nvar-1
!!                      u_1(m, l, i_dir) = 0
                      u_1(m, l, i_dir) = - u_0(m, l, i_dir)
                    enddo
                    u_1(nvar, l, i_dir) = u_0(nvar, l, i_dir)
                  enddo
                enddo

              endif

! +++++Probleme d'absence de maille voisine !

            else

              Print*
              Print*, ' Flux_visc : maille suivante absente !!', &
                      ' direction = ', dir, &
                      ' niv = ', courant%maille_niv, &
                      ' indices = ', courant%maille_i
              Stop

            endif

! +++++Calcul des composantes du flux visqueux suivant direction "dir"

            mu_ip12 = 0.5 * ( mu_1 + mu_0 )
!!            mu_ip12 = 2. * mu_1 * mu_0 / ( mu_1 + mu_0 )

            do m = 1, ndim
              u_ip12(m) = 0.5 * ( u_1(m+1, 0, dir)/u_1(1, 0, dir) + &
                                  u_0(m+1, 0, dir)/u_0(1, 0, dir) )
            enddo

            do m = 2, nvar-1         
              u_xyz(m-1, dir) = ( u_1(m, 0, dir)/u_1(1, 0, dir) - &
                                  u_0(m, 0, dir)/u_0(1, 0, dir) ) &
                                / maille%dx(dir)
            enddo
          
            ec = 0.
            do m = 2, nvar-1
              ec = ec +  u_1(m, 0, dir)**2
            enddo
            ec = 0.5 * ec / u_1(1, 0, dir)**2

            T_1 = gamma*(gamma-1.)*Mach**2 * &
                  ( u_1(nvar, 0, dir)/u_1(1, 0, dir) - ec )

            ec = 0.
            do m = 2, nvar-1
              ec = ec +  u_0(m, 0, dir)**2
            enddo
            ec = 0.5 * ec / u_0(1, 0, dir)**2

            T_0 = gamma*(gamma-1.)*Mach**2 * &
                  ( u_0(nvar, 0, dir)/u_0(1, 0, dir) - ec )

            T_xyz = ( T_1 - T_0 ) / maille%dx(dir)

! .....Directions perpendiculaires

            if( ndim > 1 ) then

              do l = i_beg(dir), i_end(dir), i_pas(dir) 

                do m = 2, nvar-1
                  u_xyz(m-1, l) = 0.25 * ( u_1(m, 1, l) / u_1(1, 1, l) - &
                                           u_1(m,-1, l) / u_1(1,-1, l) &
                                         + u_0(m, 1, l) / u_0(1, 1, l) - &
                                           u_0(m,-1, l) / u_0(1,-1, l) &
                                          ) / maille%dx(l)
                enddo            
          
              enddo

            endif          

            divV = 0.
            do m = 2, nvar-1
              divV = divV + u_xyz(m-1,m-1) 
            enddo
          
            C_visc(1) = 0.
            do m = 2, nvar-1
              C_visc(m) = - mu_ip12*( u_xyz(dir,m-1) + u_xyz(m-1,dir) - &
                               2.*uns3*divV*delta_ij(dir,m-1) ) / Reynolds
            enddo 

            C_visc(nvar) = 0.
            do l = 1, ndim
              C_visc(nvar) = C_visc(nvar) + u_ip12(l) * C_visc(l+1) 
            enddo
            C_visc(nvar) = C_visc(nvar) -  &
                  mu_ip12 * T_xyz / ( (gamma-1.)*Prandtl*Reynolds*Mach**2 )

! +++++Somme des flux pour conservation

            fvisc_droit(:) = fvisc_droit(:) + C_visc(:) / float(2**(ndim-1))

! +++++Transfert des flux pour le flux a gauche de la feuille suivante 

            local => courant%support(dir)%suivant
            local%fvisc_gch(:,dir) = C_visc(:)

          enddo

! +++++Sinon, on reconstruit au niveau de la feuille

        else

! ======== Flux visqueux a droite pour toutes les feuilles =====

          fvisc_droit(:) = 0.

! *****Maille courante

          courant => maille

! .....vicosite dynamique et conductivite thermique
          call viscosite(courant, mu_0)

! .....vitesses au centre de la maille

          do l = 1, ndim
            do m = 1, nvar
              u_0(m, 0, l) = courant%u(m)
            enddo
          enddo

! *****Directions perpendiculaires

          if( ndim > 1 ) then

            do i_dir = i_beg(dir), i_end(dir), i_pas(dir)

! +++++Point precedent dans les directions perpendiculaires

              if( associated(courant%support(i_dir)%precedent) ) then

                local => courant%support(i_dir)%precedent

                do m = 1, nvar
                  u_0(m, -1, i_dir) = local%u(m)
                enddo

! +++++Conditions aux limites (sauf periodiques)

              elseif( courant%maille_i(i_dir) == 1 ) then

! .....Conditions de symetrie

                if( caltype_deb(i_dir) == -1 ) then

                  do m = 1, nvar
                    u_0(m, -1, i_dir) = u_0(m, 0, i_dir)
                  enddo
                  u_0(i_dir+1, -1, i_dir) = - u_0(i_dir+1, 0, i_dir)

! .....Conditions de frontiere fluide ==> derivee nulle
  
                elseif( caltype_deb(i_dir) == 0 ) then

                  do m = 1, nvar
                    u_0(m, -1, i_dir) = u_0(m, 0, i_dir)
                  enddo

! .....Conditions de paroi

                elseif( caltype_deb(i_dir) == 1 ) then

                  u_0(1, -1, i_dir) = u_0(1, 0, i_dir)
                  do m = 2, nvar-1
!!                    u_0(m, -1, i_dir) = 0
                    u_0(m, -1, i_dir) = - u_0(m, 0, i_dir)
                  enddo
                  u_0(nvar, -1, i_dir) = u_0(nvar, 0, i_dir)

                endif

! +++++Probleme d'absence de maille voisine !

              else

                Print*
                Print*, ' Flux_visc : maille precedente absente !!', &
                        ' direction = ', dir,' transverse = ',i_dir, &
                        ' niv = ', courant%maille_niv, &
                          ' indices = ', courant%maille_i, &
                          ' niv = ', local%maille_niv, &
                          ' indices = ', local%maille_i

                Stop

              endif

! +++++Point suivant dans les directions perpendiculaires

              if( associated(courant%support(i_dir)%suivant) ) then

                local => courant%support(i_dir)%suivant

                do m = 1, nvar
                  u_0(m, 1, i_dir) = local%u(m)
                enddo

! +++++Conditions aux limites

              elseif( courant%maille_i(i_dir) == &
                      2**courant%maille_niv * nbre_racine(i_dir) ) then

! .....Conditions de symetrie

                if( caltype_fin(i_dir) == -1 ) then

                  do m = 1, nvar
                    u_0(m, 1, i_dir) = u_0(m, 0, i_dir)
                  enddo
                  u_0(i_dir+1, 1, i_dir) = - u_0(i_dir+1, 0, i_dir)

! .....Conditions de frontiere fluide ==> derivee nulle
  
                elseif( caltype_fin(i_dir) == 0 ) then

                  do m = 1, nvar
                   u_0(m, 1, i_dir) = u_0(m, 0, i_dir)
                  enddo

! .....Conditions de paroi

                elseif( caltype_fin(i_dir) == 1 ) then

                  u_0(1, 1, i_dir) = u_0(1, 0, i_dir)
                  do m = 2, nvar-1
!!                    u_0(m, 1, i_dir) = 0
                    u_0(m, 1, i_dir) = - u_0(m, 0, i_dir)
                  enddo
                  u_0(nvar, 1, i_dir) = u_0(nvar, 0, i_dir)

                endif

! +++++Probleme d'absence de maille voisine !

              else

                Print*
                Print*, ' Flux_visc : maille suivante absente !!', &
                        ' direction = ', dir,' transverse = ', i_dir, &
                        ' niv = ', courant%maille_niv, &
                          ' indices = ', courant%maille_i, &
                          ' niv = ', local%maille_niv, &
                          ' indices = ', local%maille_i

                Stop

              endif

            enddo

          endif

! *****Maille suivante

          if( associated(courant%support(dir)%suivant) ) then

            crt_suivant => courant%support(dir)%suivant

! .....vicosite dynamique et conductivite thermique
            call viscosite(crt_suivant, mu_1)

! .....vitesses au centre de la maille
            do l = 1, ndim
              do m = 1, nvar
                u_1(m, 0, l) = crt_suivant%u(m)
              enddo
            enddo

! *****Directions perpendiculaires

            if( ndim > 1 ) then

              do i_dir = i_beg(dir), i_end(dir), i_pas(dir)

! +++++Point precedent dans les directions perpendiculaires

                if( associated(crt_suivant%support(i_dir)%precedent) ) then

                  local => crt_suivant%support(i_dir)%precedent

                  do m = 1, nvar
                    u_1(m, -1, i_dir) = local%u(m)
                  enddo

! +++++Conditions aux limites

                elseif( crt_suivant%maille_i(i_dir) == 1 ) then

! .....Conditions de symetrie

                  if( caltype_deb(i_dir) == -1 ) then

                    do m = 1, nvar
                      u_1(m, -1, i_dir) = u_1(m, 0, i_dir)
                    enddo
                    u_1(i_dir+1, -1, i_dir) = - u_1(i_dir+1, 0, i_dir)

! .....Conditions de frontiere fluide ==> derivee nulle
  
                  elseif( caltype_deb(i_dir) == 0 ) then

                    do m = 1, nvar
                      u_1(m, -1, i_dir) = u_1(m, 0, i_dir)
                    enddo

! .....Conditions de paroi

                  elseif( caltype_deb(i_dir) == 1 ) then

                    u_1(1, -1, i_dir) = u_1(1, 0, i_dir)
                    do m = 2, nvar-1
!!                      u_1(m, -1, i_dir) = 0
                      u_1(m, -1, i_dir) = - u_1(m, 0, i_dir)
                    enddo
                    u_1(nvar, -1, i_dir) = u_1(nvar, 0, i_dir)

                  endif

! +++++Probleme d'absence de maille voisine !

                else

                  Print*
                  Print*, ' Flux_visc : maille precedente absente !!', &
                          ' direction = ',dir,' transverse = ', i_dir, &
                          ' niv = ', crt_suivant%maille_niv, &
                          ' indices = ', courant%maille_i, &
                          ' niv = ', local%maille_niv, &
                          ' indices = ', local%maille_i

                  Stop

                endif

! +++++Point suivant dans les directions perpendiculaires

                if( associated(crt_suivant%support(i_dir)%suivant) ) then

                  local => crt_suivant%support(i_dir)%suivant

                  do m = 1, nvar
                    u_1(m, 1, i_dir) = local%u(m)
                  enddo

! +++++Conditions aux limites

                elseif( crt_suivant%maille_i(i_dir) == &
                        2**crt_suivant%maille_niv * nbre_racine(i_dir) ) then

! .....Conditions de symetrie

                  if( caltype_fin(i_dir) == -1 ) then

                    do m = 1, nvar
                      u_1(m, 1, i_dir) = u_1(m, 0, i_dir)
                    enddo
                    u_1(i_dir+1, 1, i_dir) = - u_1(i_dir+1, 0, i_dir)

! .....Conditions de frontiere fluide ==> derivee nulle
  
                  elseif( caltype_fin(i_dir) == 0 ) then

                    do m = 1, nvar
                      u_1(m, 1, i_dir) = u_1(m, 0, i_dir)
                    enddo

! .....Conditions de paroi

                  elseif( caltype_fin(i_dir) == 1 ) then

                    u_1(1, 1, i_dir) = u_1(1, 0, i_dir)
                    do m = 2, nvar-1
!!                      u_1(m, 1, i_dir) = 0
                      u_1(m, 1, i_dir) = - u_1(m, 0, i_dir)
                    enddo
                    u_1(nvar, 1, i_dir) = u_1(nvar, 0, i_dir)

                  endif

! +++++Probleme d'absence de maille voisine !

                else

                  Print*
                  Print*, ' Flux_visc : maille suivante absente !!', &
                          ' direction = ',dir,' transverse = ', i_dir, &
                          ' niv = ', crt_suivant%maille_niv, &
                          ' indices = ', courant%maille_i, &
                          ' niv = ', local%maille_niv, &
                          ' indices = ', local%maille_i

                  Stop

                endif

              enddo

            endif

! +++++Conditions aux limites suivant dir

          elseif( courant%maille_i(dir) == &
                  2**courant%maille_niv * nbre_racine(dir) ) then

            mu_1 = mu_0

! .....Conditions de symetrie

            if( caltype_fin(dir) == -1 ) then

              do i_dir = 1, ndim
                do l = -1, 1
                  do m = 1, nvar
                    u_1(m, l, i_dir) = u_0(m, l, i_dir)
                  enddo
                  u_1(dir+1, l, i_dir) = - u_0(dir+1, l, i_dir)
                enddo
              enddo

! .....Conditions de frontiere fluide ==> derivee nulle
  
            elseif( caltype_fin(dir) == 0 ) then

              do i_dir = 1, ndim
                do l = -1, 1
                  do m = 1, nvar
                    u_1(m, l, i_dir) = u_0(m, l, i_dir)
                  enddo
                enddo
              enddo

! .....Conditions de paroi

            elseif( caltype_fin(dir) == 1 ) then

              do i_dir = 1, ndim
                do l = -1, 1
                  u_1(1, l, i_dir) = u_0(1, l, i_dir)
                  do m = 2, nvar-1
!!                    u_1(m, l, i_dir) = 0
                    u_1(m, l, i_dir) = - u_0(m, l, i_dir)
                  enddo
                  u_1(nvar, l, i_dir) = u_0(nvar, l, i_dir)
                enddo
              enddo

            endif

! +++++Probleme d'absence de maille voisine !

          else

            Print*
            Print*, ' Flux_visc : maille suivante absente !!', &
                    ' direction = ', dir, &
                    ' niv = ', maille%maille_niv, &
                   ' indices = ', maille%maille_i
            Stop

          endif

! +++++Calcul des composantes du flux visqueux suivant X

          mu_ip12 = 0.5 * ( mu_1 + mu_0 )
!!          mu_ip12 = 2. * mu_1 * mu_0 / ( mu_1 + mu_0 )

          do m = 1, ndim
            u_ip12(m) = 0.5 * ( u_1(m+1, 0, dir) / u_1(1, 0, dir) + &
                                u_0(m+1, 0, dir) / u_0(1, 0, dir) )
          enddo

          do m = 2, nvar-1         
             u_xyz(m-1, dir) = ( u_1(m, 0, dir) / u_1(1, 0, dir) - &
                                 u_0(m, 0, dir) / u_0(1, 0, dir) ) &
                               / maille%dx(dir)
          enddo
          
          ec = 0.
          do m = 2, nvar-1
            ec = ec +  u_1(m, 0, dir)**2
          enddo
          ec = 0.5 * ec / u_1(1, 0, dir)**2

          T_1 = gamma*(gamma-1.)*Mach**2 * &
                ( u_1(nvar, 0, dir)/u_1(1, 0, dir) - ec )

          ec = 0.
          do m = 2, nvar-1
            ec = ec +  u_0(m, 0, dir)**2
          enddo
          ec = 0.5 * ec / u_0(1, 0, dir)**2

          T_0 = gamma*(gamma-1.)*Mach**2 * &
                ( u_0(nvar, 0, dir)/u_0(1, 0, dir) - ec )

          T_xyz = ( T_1 - T_0 ) / maille%dx(dir)

! .....Directions perpendiculaires

          if( ndim > 1 ) then

            do l = i_beg(dir), i_end(dir), i_pas(dir)

              do m = 2, nvar-1
                u_xyz(m-1, l) = 0.25 * ( u_1(m, 1, l) / u_1(1, 1, l) - &
                                         u_1(m,-1, l) / u_1(1,-1, l) &
                                       + u_0(m, 1, l) / u_0(1, 1, l) - &
                                         u_0(m,-1, l) / u_0(1,-1, l) &
                                        ) / maille%dx(l)
              enddo            
          
            enddo

          endif          

          divV = 0.
          do m = 2, nvar-1
            divV = divV + u_xyz(m-1,m-1) 
          enddo
          
          C_visc(1) = 0.
          do m = 2, nvar-1
            C_visc(m) = - mu_ip12*( u_xyz(dir,m-1) + u_xyz(m-1,dir) - &
                             2.*uns3*divV*delta_ij(dir,m-1) ) / Reynolds
          enddo 

          C_visc(nvar) = 0.
          do l = 1, ndim
            C_visc(nvar) = C_visc(nvar) + u_ip12(l) * C_visc(l+1) 
          enddo
          C_visc(nvar) = C_visc(nvar) - &
              mu_ip12 * T_xyz / ( (gamma-1.)*Prandtl*Reynolds*Mach**2 )

! +++++Somme des flux pour conservation

          fvisc_droit(:) = C_visc(:)

! +++++Transfert du flux pour le flux a gauche de la feuille suivante,
!      si elle existe 

          if( associated(maille%support(dir)%suivant) ) &
            maille%support(dir)%suivant%fvisc_gch(:,dir) = C_visc(:) 

        endif

! .....Flux Visqueux a droite de la maille courante

        maille%fvisc_dt(:,dir) = fvisc_droit(:)

      Enddo

101   format(' Niv > : ',3(1x,i6))
102   format(4(1x,1pe15.8))
103   format(' Niv = : ',3(1x,i6))

      end subroutine flux_visc
