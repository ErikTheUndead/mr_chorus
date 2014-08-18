      subroutine flux_osmp7_conv(dir)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: dir

      LOGICAL :: flag

      INTEGER :: i_sym, l, l_dt, leaf, m
      INTEGER :: s_deb, s_fin
      INTEGER :: ijk_gch
      INTEGER, DIMENSION(4,3) :: ijk_dt

      REAL(kind=8), DIMENSION(nvar, -nordre+1:nordre) :: u_n
      REAL(kind=8), DIMENSION(nvar) :: C_mp, flux_conv
      REAL(kind=8), DIMENSION(nvar) :: fmp_droit,f_droit

      REAL(kind=8) :: rho_ec, sigma
      REAL(kind=8) :: rho_wall, rhoE_wall

      TYPE(structure_maille), pointer :: maille
      TYPE(structure_maille), pointer :: local, courant

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

! ======== Correction OSMP7 : Condition a la limite amont du domaine =====

        if(maille%maille_i(dir) == 1 .and. caltype_deb(dir) /= 2) then

          local => maille

! .....Rapport : sigma/2 = dt/2 / dx
          sigma = dts2 / local%dx(dir)

! +++++Pour toutes les mailles utiles

          s_deb = 1
          s_fin = 1

          do l = s_deb, nordre

            do m = 1, nvar
              u_n(m,l) = local%u(m)
            enddo

! .....on passe a la maille suivante
            s_fin = l

!!CT0110            if( .not.associated(local%support(dir)%suivant) ) exit
            local => local%support(dir)%suivant

          enddo

! +++++Reconstruction en amont du domaine

! .....Conditions de symetrie

          if( caltype_deb(dir) == -1 ) then

            do l = -nordre+1, s_deb-1

              i_sym = 2*s_deb-(l+1)

              do m = 1, nvar
                u_n(m,l) = u_n(m,i_sym)
              enddo
              u_n(dir+1,l) = -u_n(dir+1,i_sym)

            enddo

! .....Conditions de frontiere fluide ==> derivee nulle // Neumann

          elseif( caltype_deb(dir) == 0 ) then

            do l = -nordre+1, s_deb-1

              i_sym = 2*s_deb-(l+1)

               u_n(1,l) = u_n(1, s_deb)
               do m = 2, nvar-1
                  u_n(m,l) = v_right(m-1)*u_n(1,l)
               end do
               u_n(nvar,l) = u_n(nvar, s_deb)

!               do m = 1, nvar
! !!                u_n(m,l) = u_n(m,i_sym)
!                 u_n(m,l) = u_n(m,s_deb)
! !!                u_n(m,l) = u_n(m,s_deb) - &
! !!                           float(s_deb-l)*(u_n(m,s_deb+1)-u_n(m,s_deb))
!               enddo

            enddo

! .....Conditions de paroi

          elseif( caltype_deb(dir) == 1 ) then

            if( cle_twall ) then
! .....paroi adiabatique ==> gradient densite et presssion  = 0
              rho_wall = (9.*u_n(1,s_deb) - u_n(1,s_deb+1))/8.

              rho_ec = 0.
              do m = 2, nvar-1
                rho_ec = rho_ec + u_n(m,s_deb)**2
              enddo

              rhoE_wall = 9.*(u_n(nvar,s_deb) - 0.5*rho_ec/u_n(1,s_deb))

              rho_ec = 0.
              do m = 2, nvar-1
                rho_ec = rho_ec + u_n(m,s_deb+1)**2
              enddo

              rhoE_wall = (rhoE_wall  - &
                           (u_n(nvar,s_deb+1) - 0.5*rho_ec/u_n(1,s_deb+1)) )/8.

              do l = -nordre+1, s_deb-1
                i_sym = 2*s_deb-(l+1)

                u_n(1, l) = rho_wall
                do m = 2, nvar-1
                    u_n(m, l) = v_left(m-1)*u_n(1,l)
!                   u_n(m, l) = - u_n(m, i_sym)
                enddo
                u_n(nvar, l) = rhoE_wall
              enddo

            else
! .....paroi  thermostatee

              do l = -nordre+1, s_deb-1

                i_sym = 2*s_deb-(l+1)

!!                u_n(1, l) = u_n(1, s_deb)
                u_n(1, l) = u_n(1, i_sym)
                do m = 2, nvar-1
                  u_n(m, l) = - u_n(m, i_sym)
!                   u_n(m,l) = v_left(m-1)*u_n(1,l)
                enddo
!!                u_n(nvar, l) = u_n(nvar, s_deb)
                u_n(nvar, l) = u_n(nvar, i_sym)
              enddo
            endif

          endif

! +++++Calcul de la correction OSMP7 a la frontiere amont du domaine : Cj+1/2

          call osmp7_conv(dir, sigma, u_n, C_mp, flux_conv)

! .....Correction TVD-MP a gauche de la maille courante, Cj+1/2 et Fj+1/2

          maille%fmp_gch(:) = C_mp(:) ! Cj+1/2
          maille%f_gch(:) = flux_conv(:) ! Fj+1/2

        endif

! ============ Correction OSMP7 a droite pour toutes les feuilles ========

        flag = .false.

! .....Test pour reconstruire niveau feuille ou niveau sup ?

        if( associated(maille%support(dir)%suivant) ) then

          if( associated(maille%support(dir)%suivant%fils(1)%ptr) ) then

            if(maille%support(dir)%suivant%fils(1)%ptr%feuille) &
              flag = .true.

          endif
        endif

! *****Si la feuille a des fils et que
!      la feuille suivante a des fils non-fictifs

        if( flag .and. associated(maille%fils(1)%ptr) ) then

          fmp_droit(:) = 0.
          f_droit(:) = 0.

! +++++La correction est evaluee au niveau superieur (Fils de la feuille)

          do l_dt = 1, 2**(ndim-1)

            local => maille%fils( ijk_dt(l_dt,dir) )%ptr
            courant => maille%fils( ijk_dt(l_dt,dir) )%ptr

! .....Rapport : sigma/2 = dt/2 / dx
            sigma = dts2 / local%dx(dir)

! *****initialisation pour toutes les mailles utiles

            do l = 0, -nordre+1, -1
              s_deb = l
              if( .not.associated(local%support(dir)%precedent) &
                  .or. l == -nordre+1 ) exit
              local => local%support(dir)%precedent
            enddo

! .....Pour toutes les mailles utiles

            s_fin = s_deb
            do l = s_deb, nordre

              do m = 1, nvar
                u_n(m,l) = local%u(m)
              enddo

! .....on passe a la maille suivante
              s_fin = l

              if( .not.associated(local%support(dir)%suivant) ) exit
              local => local%support(dir)%suivant

            enddo

! *****Conditions a la limite amont

            if( s_deb > -nordre+1 ) then

! .....Conditions de symetrie

              if( caltype_deb(dir) == -1 ) then

                do l = -nordre+1, s_deb-1

                  i_sym = 2*s_deb-(l+1)

                  do m = 1, nvar
                    u_n(m,l) = u_n(m,i_sym)
                  enddo
                  u_n(dir+1,l) = -u_n(dir+1,i_sym)

                enddo

! .....Conditions de frontiere fluide ==> derivee nulle

              elseif( caltype_deb(dir) == 0 ) then

                do l = -nordre+1, s_deb-1

!                   i_sym = 2*s_deb-(l+1)

                     u_n(1,l) = u_n(1, s_deb)
                     do m = 2, nvar-1
                        u_n(m,l) = v_right(m-1)*u_n(1,l)
                     end do
                     u_n(nvar,l) = u_n(nvar, s_deb)


!                   do m = 1, nvar
! !!                    u_n(m,l) = u_n(m,i_sym)
!                     u_n(m,l) = u_n(m,s_deb)
! !!                    u_n(m,l) = u_n(m,s_deb) - &
! !!                               float(s_deb-l)*(u_n(m,s_deb+1)-u_n(m,s_deb))
!                   enddo

                enddo

! .....Conditions de paroi

              elseif( caltype_deb(dir) == 1 ) then

                if( cle_twall ) then
! .....paroi adiabatique ==> gradient densite et presssion  = 0
                  rho_wall = (9.*u_n(1,s_deb) - u_n(1,s_deb+1))/8.

                  rho_ec = 0.
                  do m = 2, nvar-1
                    rho_ec = rho_ec + u_n(m,s_deb)**2
                  enddo

                  rhoE_wall = 9.*(u_n(nvar,s_deb) - 0.5*rho_ec/u_n(1,s_deb))

                  rho_ec = 0.
                  do m = 2, nvar-1
                    rho_ec = rho_ec + u_n(m,s_deb+1)**2
                  enddo

                  rhoE_wall = (rhoE_wall  - &
                           (u_n(nvar,s_deb+1) - 0.5*rho_ec/u_n(1,s_deb+1)) )/8.

                  do l = -nordre+1, s_deb-1
                    i_sym = 2*s_deb-(l+1)

                    u_n(1, l) = rho_wall
                    do m = 2, nvar-1
                     u_n(m, l) = v_left(m-1)*u_n(1,l)
!                       u_n(m,l) = - u_n(m, i_sym)
                    enddo
                    u_n(nvar, l) = rhoE_wall
                  enddo

                else
! .....paroi  thermostatee

                  do l = -nordre+1, s_deb-1

                    i_sym = 2*s_deb-(l+1)

!!                    u_n(1,l) = u_n(1,s_deb)
                    u_n(1,l) = u_n(1, i_sym)
                    do m = 2, nvar-1
                      u_n(m,l) = - u_n(m, i_sym)
!                         u_n(m,l) = v_left(m-1)*u_n(1,l)
                    enddo
!!                    u_n(nvar,l) = u_n(nvar,s_deb)
                    u_n(nvar,l) = u_n(nvar, i_sym)
                  enddo

                endif
              endif

            endif

! *****Conditions a la limite aval

            if( s_fin < nordre ) then

! .....Conditions de symetrie

              if( caltype_fin(dir) == -1 ) then

                do l = s_fin+1, nordre

                  i_sym = 2*s_fin-(l-1)

                  do m = 1, nvar
                    u_n(m,l) = u_n(m,i_sym)
                  enddo
                  u_n(dir+1,l) = -u_n(dir+1,i_sym)

                enddo

! .....Conditions de frontiere fluide ==> derivee nulle

              elseif( caltype_fin(dir) == 0 ) then

                do l = s_fin+1, nordre

!                   i_sym = 2*s_fin-(l-1)

                     u_n(1,l) = u_n(1, s_fin)
!                      u_n(2,l) = v_right(1)*u_n(1,l)
!                      u_n(3,l) = u_n(3,s_fin)
                     do m = 2, nvar-1
                        u_n(m,l) = v_right(m-1)*u_n(1,l)
                     end do
                     u_n(nvar,l) = u_n(nvar,s_fin)
!                   do m = 1, nvar
! !!                    u_n(m,l) = u_n(m,i_sym)
!                     u_n(m,l) = u_n(m,s_fin)
! !!                    u_n(m,l) = u_n(m,s_fin) + &
! !!                               float(l-s_fin)*(u_n(m,s_fin)-u_n(m,s_fin-1))
!                   enddo

                enddo

! .....Conditions de paroi

              elseif( caltype_fin(dir) == 1 ) then

                if( cle_twall ) then
! .....paroi adiabatique ==> gradient densite et presssion  = 0
                  rho_wall = (9.*u_n(1,s_fin) - u_n(1,s_fin-1))/8.

                  rho_ec = 0.
                  do m = 2, nvar-1
                    rho_ec = rho_ec + u_n(m,s_fin)**2
                  enddo

                  rhoE_wall = 9.*(u_n(nvar,s_fin) - 0.5*rho_ec/u_n(1,s_fin))

                  rho_ec = 0.
                  do m = 2, nvar-1
                    rho_ec = rho_ec + u_n(m,s_fin-1)**2
                  enddo

                  rhoE_wall = (rhoE_wall  - &
                           (u_n(nvar,s_fin-1) - 0.5*rho_ec/u_n(1,s_fin-1)) )/8.

                  do l = s_fin+1, nordre
                    i_sym = 2*s_fin-(l-1)

                    u_n(1, l) = rho_wall
                    do m = 2, nvar-1
                     u_n(m, l) = v_left(m-1)*u_n(1,l)
!                       u_n(m,l) = - u_n(m, i_sym)
                    enddo
                    u_n(nvar, l) = rhoE_wall
                  enddo

                else
! .....paroi  thermostatee

                  do l = s_fin+1, nordre

                    i_sym = 2*s_fin-(l-1)

!!                    u_n(1,l) = u_n(1,s_fin)
                    u_n(1,l) = u_n(1, i_sym)
                    do m = 2, nvar-1
                      u_n(m,l) = - u_n(m, i_sym)
!                         u_n(m, l) = v_left(m-1)*u_n(1,l)
                    enddo
!!                    u_n(nvar,l) = u_n(nvar,s_fin)
                    u_n(nvar,l) = u_n(nvar, i_sym)
                  enddo
                endif

              endif

            endif

! +++++Calcul de la correction OSMP7 a droite

            call osmp7_conv(dir, sigma, u_n, C_mp, flux_conv)

! .....Somme des flux pour conservation

            fmp_droit(:) = fmp_droit(:) + &
                            C_mp(:) / float(2**(ndim-1))

            f_droit(:) = f_droit(:) + &
                            flux_conv(:) / float(2**(ndim-1))

! +++++Transfert des flux pour le flux a gauche de la feuille suivante

! .....indice de la maille suivante
            ijk_gch = ijk_dt( l_dt,dir ) - 2**(dir-1)

            maille%support(dir)%suivant%fils(ijk_gch)%ptr%fmp_gch(:)= C_mp(:)
            maille%support(dir)%suivant%fils(ijk_gch)%ptr%f_gch(:)= flux_conv(:)

          enddo

! .....Correction TVD-MP a droite pour la maille courante

          maille%fmp_dt(:) = fmp_droit(:)
          maille%f_dt(:) = f_droit(:)

        else

! +++++Sinon, on reconstruit au niveau de la feuille

          local => maille
          courant => maille

! .....Rapport : sigma/2 = dt/2 / dx
          sigma = dts2 / local%dx(dir)

! .....initialisation du support

          do l = 0, -nordre+1, -1
            s_deb = l
            if( .not.associated(local%support(dir)%precedent) &
                .or. l == -nordre+1 ) exit
            local => local%support(dir)%precedent
          enddo

! *****Pour toutes les mailles utiles

          s_fin = s_deb
          do l = s_deb, nordre

            do m = 1, nvar
              u_n(m,l) = local%u(m)
            enddo

! .....on passe a la maille suivante
            s_fin = l

            if( .not.associated(local%support(dir)%suivant) ) exit
            local => local%support(dir)%suivant

          enddo

! *****Conditions a la limite amont

          if( s_deb > -nordre+1 ) then

! .....Conditions de symetrie

            if( caltype_deb(dir) == -1 ) then

              do l = -nordre+1, s_deb-1

                i_sym = 2*s_deb-(l+1)

                do m = 1, nvar
                  u_n(m,l) = u_n(m,i_sym)
                enddo
                u_n(dir+1,l) = -u_n(dir+1,i_sym)

              enddo

! .....Conditions de frontiere fluide ==> derivee nulle

            elseif( caltype_deb(dir) == 0 ) then

              do l = -nordre+1, s_deb-1

!                 i_sym = 2*s_deb-(l+1)

                  u_n(1,l) = u_n(1, s_deb)
                  do m = 2, nvar-1
                     u_n(m,l) = v_right(m-1)*u_n(1,l)
                  end do
                  u_n(nvar,l) = u_n(nvar, s_deb)

!                 do m = 1, nvar
! !!                  u_n(m,l) = u_n(m,i_sym)
!                   u_n(m,l) = u_n(m, s_deb)
! !!                  u_n(m,l) = u_n(m,s_deb) - &
! !!                             float(s_deb-l)*(u_n(m,s_deb+1)-u_n(m,s_deb))
!                 enddo

              enddo

! .....Conditions de paroi

            elseif( caltype_deb(dir) == 1 ) then

              if( cle_twall ) then
! .....paroi adiabatique ==> gradient densite et presssion  = 0
                rho_wall = (9.*u_n(1,s_deb) - u_n(1,s_deb+1))/8.

                rho_ec = 0.
                do m = 2, nvar-1
                  rho_ec = rho_ec + u_n(m,s_deb)**2
                enddo

                rhoE_wall = 9.*(u_n(nvar,s_deb) - 0.5*rho_ec/u_n(1,s_deb))

                rho_ec = 0.
                do m = 2, nvar-1
                  rho_ec = rho_ec + u_n(m,s_deb+1)**2
                enddo

                rhoE_wall = (rhoE_wall  - &
                         (u_n(nvar,s_deb+1) - 0.5*rho_ec/u_n(1,s_deb+1)) )/8.

                do l = -nordre+1, s_deb-1
                  i_sym = 2*s_deb-(l+1)

                  u_n(1, l) = rho_wall
                  do m = 2, nvar-1
                   u_n(m, l) = v_left(m-1)*u_n(1,l)
!                     u_n(m,l) = - u_n(m, i_sym)
                  enddo
                  u_n(nvar, l) = rhoE_wall
                enddo

              else
! .....paroi  thermostatee

                do l = -nordre+1, s_deb-1

                  i_sym = 2*s_deb-(l+1)

!!                  u_n(1,l) = u_n(1,s_deb)
                  u_n(1,l) = u_n(1, i_sym)
                  do m = 2, nvar-1
                    u_n(m,l) = - u_n(m, i_sym)
!                      u_n(m,l) = v_left(m-1)*u_n(1,l)
                  enddo
!!                  u_n(nvar,l) = u_n(nvar,s_deb)
                  u_n(nvar,l) = u_n(nvar, i_sym)
                enddo
              endif
            endif

          endif

! *****Conditions a la limite aval

          if( s_fin < nordre ) then

! .....Conditions de symetrie

            if( caltype_fin(dir) == -1 ) then

              do l = s_fin+1, nordre

                i_sym = 2*s_fin-(l-1)

                do m = 1, nvar
                  u_n(m,l) = u_n(m,i_sym)
                enddo
                u_n(dir+1,l) = -u_n(dir+1,i_sym)

              enddo

! .....Conditions de frontiere fluide ==> derivee nulle

            elseif( caltype_fin(dir) == 0 ) then

              do l = s_fin+1, nordre

!                 i_sym = 2*s_fin-(l-1)

                  u_n(1,l) = u_n(1,s_fin)
!                   u_n(2,l) = v_right(1)*u_n(1,l)
!                   u_n(3,l) = u_n(3,s_fin)
                  do m = 2, nvar-1
                     u_n(m,l) = v_right(m-1)*u_n(1,l)
                  end do
                  u_n(nvar,l) = u_n(nvar, s_fin)
!                 do m = 1, nvar
! !!                  u_n(m,l) = u_n(m,i_sym)
!                   u_n(m,l) = u_n(m, s_fin)
! !!                  u_n(m,l) = u_n(m,s_fin) + &
! !!                             float(l-s_fin)*(u_n(m,s_fin)-u_n(m,s_fin-1))
!                 enddo

              enddo

! .....Conditions de paroi

            elseif( caltype_fin(dir) == 1 ) then

              if( cle_twall ) then
! .....paroi adiabatique ==> gradient densite et presssion  = 0
                rho_wall = (9.*u_n(1,s_fin) - u_n(1,s_fin-1))/8.

                rho_ec = 0.
                do m = 2, nvar-1
                  rho_ec = rho_ec + u_n(m,s_fin)**2
                enddo

                rhoE_wall = 9.*(u_n(nvar,s_fin) - 0.5*rho_ec/u_n(1,s_fin))

                rho_ec = 0.
                do m = 2, nvar-1
                  rho_ec = rho_ec + u_n(m,s_fin-1)**2
                enddo

                rhoE_wall = (rhoE_wall  - &
                         (u_n(nvar,s_fin-1) - 0.5*rho_ec/u_n(1,s_fin-1)) )/8.

                do l = s_fin+1, nordre
                  i_sym = 2*s_fin-(l-1)

                  u_n(1, l) = rho_wall
                  do m = 2, nvar-1
                   u_n(m, l) = v_left(m-1)*u_n(1,l)
!                     u_n(m,l) = - u_n(m, i_sym)
                  enddo
                  u_n(nvar, l) = rhoE_wall
                enddo

              else
! .....paroi  thermostatee

                do l = s_fin+1, nordre

                  i_sym = 2*s_fin-(l-1)

!!                  u_n(1,l) = u_n(1,s_fin)
                  u_n(1,l) = u_n(1, i_sym)
                  do m = 2, nvar-1
                    u_n(m,l) = - u_n(m, i_sym)
!                      u_n(m,l) = v_left(m-1)*u_n(1,l)
                  enddo
!!                  u_n(nvar,l) = u_n(nvar,s_fin)
                  u_n(nvar,l) = u_n(nvar, i_sym)
                enddo
              endif

            endif

          endif

! +++++Calcul de la correction OSMP7 a droite

          call osmp7_conv(dir, sigma, u_n, C_mp, flux_conv)

! .....Correction TVD-MP a droite pour la maille courante

          maille%fmp_dt(:) = C_mp(:)
          maille%f_dt(:) = flux_conv(:)

! +++++Transfert des flux pour le flux a gauche de la feuille suivante

          if( associated(maille%support(dir)%suivant) ) then
            maille%support(dir)%suivant%fmp_gch(:) = C_mp(:)
            maille%support(dir)%suivant%f_gch(:) = flux_conv(:)
          end if

        endif

      Enddo

      end subroutine flux_osmp7_conv
