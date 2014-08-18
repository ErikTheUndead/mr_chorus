      recursive subroutine feuille_fictive_visc(maille)

      use mod_common
      use mod_structure
      use mod_recherche

      implicit none

      TYPE(structure_maille), pointer :: maille

      LOGICAL :: flag, reelle

      INTEGER :: i_dir, l, m, s_deb, st
      INTEGER, DIMENSION(3) :: i_beg, i_end, i_pas

      TYPE(structure_maille), pointer :: courant, stencil, local, pere

! +++++Initialisation des indices

      i_beg(1) = 2
      i_end(1) = ndim
      i_pas(1) = 1

      i_beg(2) = 1
      i_end(2) = ndim
      i_pas(2) = 2

      i_beg(3) = 1
      i_end(3) = 2
      i_pas(3) = 1

! +++++On marque les feuilles reelles (non-fictives)

      flag = .false.

      if( .not.associated( maille%fils(1)%ptr) ) then
        if( maille%feuille ) flag = .true.
      else
        if( .not.maille%fils(1)%ptr%feuille ) flag = .true.
      endif


! +++++On parcours les feuilles reelles (non-fictives)

      If( flag ) then

! +++++Dans toutes les directions, on teste les mailles voisines

        do l = 1, ndim

! *****Maille courante

          courant => maille

! .....On balaye les deux autres directions

          stencil => courant

          do i_dir = i_beg(l), i_end(l), i_pas(l)

! .....maille precedente dans la direction normale i_dir 

            if( .not.associated(stencil%support(i_dir)%precedent) ) then

! =====Creation a partir du pere : si pas sur un bord du domaine

              if( stencil%maille_i(i_dir) /= 1 ) then

                pere => stencil%pere%support(i_dir)%precedent

! .....Creation des fils
                reelle = .false.
                call creation(pere, reelle)

! .....comptabilisation des feuilles fictives
                compteur_fictive = compteur_fictive + 2**ndim

! .....Graduation
                call graduation_local(pere, reelle)

! .....Valorisation des fils
                call decoding_moy(pere)

              endif

! .....Si la feuille existe ==> valoriser si elle est fictive
            elseif( .not.stencil%support(i_dir)%precedent%feuille ) then

              local => stencil%support(i_dir)%precedent

! .....Valorisation des fils si feuille fictive
              if( .not.local%feuille ) then
                  
                pere => local%pere

                call decoding_moy(pere)

              endif
            endif

! .....maille suivante dans la direction normale i_dir

            if( .not.associated( stencil%support(i_dir)%suivant) ) then

! =====Creation a partir du pere : si pas sur un bord du domaine

              if( stencil%maille_i(i_dir) /= &
                  2**stencil%maille_niv * nbre_racine(i_dir) ) then

                pere => stencil%pere%support(i_dir)%suivant

! .....Creation des fils
                reelle = .false.
                call creation(pere, reelle)

! .....comptabilisation des feuilles fictives
                compteur_fictive = compteur_fictive + 2**ndim

! .....Graduation
                call graduation_local(pere, reelle)

! .....Valorisation des fils
                call decoding_moy(pere)

              endif

! .....Si la feuille existe ==> valoriser si elle est fictive
            elseif( .not.stencil%support(i_dir)%suivant%feuille ) then

              local => stencil%support(i_dir)%suivant

! .....Valorisation des fils si feuille fictive
              if( .not.local%feuille ) then
                  
                pere => local%pere

                call decoding_moy(pere)

              endif

            endif

          enddo

! *****Maille precedente dans la direction l

! +++++Si la maille precedente dans la direction l n'existe pas ==> la creer
          if( .not.associated(courant%support(l)%precedent) ) then

! .....Creation a partir du pere : si pas sur un bord du domaine
            if( courant%maille_i(l) /= 1 ) then

              pere => courant%pere%support(l)%precedent

! .....Creation des fils
              reelle = .false.
              call creation(pere, reelle)

! .....comptabilisation des feuilles fictives
              compteur_fictive = compteur_fictive + 2**ndim

! .....Graduation
              call graduation_local(pere, reelle)

! .....Valorisation des fils
              call decoding_moy(pere)

            endif
          endif

! +++++La feuille precedente dans la direction l existe

          if( associated(courant%support(l)%precedent) ) then
  
            stencil => courant%support(l)%precedent

! =====Si la feuille precedente est fictive ==> la valoriser
            if( .not.stencil%feuille ) then

                pere => stencil%pere
                call decoding_moy(pere)

            endif

! .....On balaye les deux autres directions

            do i_dir = i_beg(l), i_end(l), i_pas(l)

! .....maille precedente dans la direction normale i_dir

              if( .not.associated( stencil%support(i_dir)%precedent) ) then

                if( stencil%maille_i(i_dir) /= 1 ) then

! =====Creation a partir du pere : si pas sur un bord du domaine

                  pere => stencil%pere%support(i_dir)%precedent

! .....Creation des fils
                  reelle = .false.
                  call creation(pere, reelle)

! .....comptabilisation des feuilles fictives
                  compteur_fictive = compteur_fictive + 2**ndim

! .....Graduation
                  call graduation_local(pere, reelle)

! .....Valorisation des fils
                  call decoding_moy(pere)

                endif

! .....Si la feuille existe ==> valoriser si elle est fictive
              elseif( .not.stencil%support(i_dir)%precedent%feuille ) then

                local => stencil%support(i_dir)%precedent

! .....Valorisation des fils si feuille fictive
                if( .not.local%feuille ) then
                  
                  pere => local%pere

                  call decoding_moy(pere)

                endif
              endif

! .....maille suivante dans la direction normale i_dir

              if( .not.associated( stencil%support(i_dir)%suivant) ) then

                if( stencil%maille_i(i_dir) /= &
                    2**stencil%maille_niv * nbre_racine(i_dir) ) then

! =====Creation a partir du pere : si pas sur un bord du domaine

                  pere => stencil%pere%support(i_dir)%suivant

! .....Creation des fils
                  reelle = .false.
                  call creation(pere, reelle)

! .....comptabilisation des feuilles fictives
                  compteur_fictive = compteur_fictive + 2**ndim

! .....Graduation
                  call graduation_local(pere, reelle)

! .....Valorisation des fils
                  call decoding_moy(pere)

                endif

! .....Si la feuille existe ==> valoriser si elle est fictive
              elseif( .not.stencil%support(i_dir)%suivant%feuille ) then

                local => stencil%support(i_dir)%suivant

! .....Valorisation des fils si feuille fictive
                if( .not.local%feuille ) then
                  
                  pere => local%pere

                  call decoding_moy(pere)

                endif
              endif

            enddo

          endif

! *****Maille suivante suivant la direction l

! +++++Si la maille suivante dans la direction l  n'existe pas ==> la creer
          if( .not.associated(courant%support(l)%suivant) ) then

! .....Creation a partir du pere : si pas sur un bord du domaine
            if( courant%maille_i(l) /= &
                    2**courant%maille_niv * nbre_racine(l) ) then

              pere => courant%pere%support(l)%suivant

! .....Creation des fils
              reelle = .false.
              call creation(pere, reelle)

! .....comptabilisation des feuilles fictives
              compteur_fictive = compteur_fictive + 2**ndim

! .....Graduation
              call graduation_local(pere, reelle)

! .....Valorisation des fils
              call decoding_moy(pere)

            endif
          endif

! +++++La feuille suivante dans la direction l existe

          if( associated(courant%support(l)%suivant) ) then
  
            stencil => courant%support(l)%suivant

! =====Si la feuille suivante est fictive ==> la valoriser
            if( .not.stencil%feuille ) then

                pere => stencil%pere
                call decoding_moy(pere)

            endif

! .....On balaye les deux autres directions

            do i_dir = i_beg(l), i_end(l), i_pas(l)

! .....maille precedente dans la direction normale i_dir

              if( .not.associated( stencil%support(i_dir)%precedent) ) then

                if( stencil%maille_i(i_dir) /= 1 ) then

! =====Creation a partir du pere : si pas sur un bord du domaine

                  pere => stencil%pere%support(i_dir)%precedent

! .....Creation des fils
                  reelle = .false.
                  call creation(pere, reelle)

! .....comptabilisation des feuilles fictives
                  compteur_fictive = compteur_fictive + 2**ndim

! .....Graduation
                  call graduation_local(pere, reelle)

! .....Valorisation des fils
                  call decoding_moy(pere)

                endif

! .....Si la feuille existe ==> valoriser si elle est fictive
              elseif( .not.stencil%support(i_dir)%precedent%feuille ) then

                local => stencil%support(i_dir)%precedent

! .....Valorisation des fils si feuille fictive
                if( .not.local%feuille ) then
                  
                  pere => local%pere

                  call decoding_moy(pere)

                endif
              endif

! .....maille suivante dans la direction normale i_dir

              if( .not.associated( stencil%support(i_dir)%suivant) ) then

                if( stencil%maille_i(i_dir) /= &
                    2**stencil%maille_niv * nbre_racine(i_dir) ) then

! =====Creation a partir du pere : si pas sur un bord du domaine

                  pere => stencil%pere%support(i_dir)%suivant

! .....Creation des fils
                  reelle = .false.
                  call creation(pere, reelle)

! .....comptabilisation des feuilles fictives
                  compteur_fictive = compteur_fictive + 2**ndim

! .....Graduation
                  call graduation_local(pere, reelle)

! .....Valorisation des fils
                  call decoding_moy(pere)

                endif

! .....Si la feuille existe ==> valoriser si elle est fictive
              elseif( .not.stencil%support(i_dir)%suivant%feuille ) then

                local => stencil%support(i_dir)%suivant

! .....Valorisation des fils si feuille fictive
                if( .not.local%feuille ) then
                  
                  pere => local%pere

                  call decoding_moy(pere)

                endif
              endif

            enddo

          endif

        enddo

      Else

! .....On descend vers les fils

        if( associated(maille%fils(1)%ptr) ) then
          do l = 1, 2**ndim
            call feuille_fictive_visc(maille%fils(l)%ptr)
          enddo
        endif

      Endif

      end subroutine feuille_fictive_visc
