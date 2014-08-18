      recursive subroutine feuille_fictive(maille)

      use mod_common
      use mod_structure
      use mod_recherche
      use mod_fonction

      implicit none

      TYPE(structure_maille), pointer :: maille

      LOGICAL :: flag, reelle

      INTEGER :: l, m, s_deb, st

      TYPE(structure_maille), pointer :: stencil, courant

! +++++On marque les feuilles reelles (non-fictives)

      flag = .false.

      if( .not.associated( maille%fils(1)%ptr) ) then
        if( maille%feuille ) flag = .true.
      else
        if( .not.maille%fils(1)%ptr%feuille ) flag = .true.
      endif

! +++++on parcours les feuilles reelles (non-fictives)

      If( flag ) then

! .....On teste le niveau de la maille suivante dans toutes les directions

        do l = 1, ndim

! *****Si la feuille suivante est au meme niveau ou au niveau inferieur 
!      ==> reconstruction au niveau du pere de la feuille

          courant => maille%pere

! *****Si la feuille courante a une maille voisine suivante qui a des fils
!      ==> reconstruction au niveau superieur a partir de la feuille

          if( associated(maille%support(l)%suivant) ) then

            stencil => maille%support(l)%suivant

            if( associated(stencil%fils(1)%ptr) ) then 

              if( stencil%fils(1)%ptr%feuille ) courant => maille

            endif

          endif

! +++++Reconstruction 

! .....on remonte au point le plus a gauche du stencil

          stencil => courant
          s_deb = 0

          do st = -1, -(nordre+1)/2+1, -1

            if( associated(stencil%support(l)%precedent) ) then

              stencil => stencil%support(l)%precedent
              s_deb = st

            endif

          enddo

! .....On balaye tout le stencil

          do st = s_deb, (nordre+1)/2

! .....Si non existant, on cree et valorise 
!      le point du stencil

            if(.not.associated(stencil%fils(1)%ptr)) then

! .....Creation des fils
              reelle = .false.
              call creation(stencil, reelle)

! .....comptabilisation des feuilles fictives
              compteur_fictive = compteur_fictive + 2**ndim

! .....Graduation
              call graduation_local(stencil, reelle)

! .....Valorisation des fils
              call decoding_moy(stencil)

! .....Sinon, on valorise le point du stencil

            else

! .....Valorisation des fils si fictifs
              if( .not.stencil%fils(1)%ptr%feuille ) &
                call decoding_moy(stencil)

            endif

            if( .not.associated(stencil%support(l)%suivant) ) goto 10
            stencil => stencil%support(l)%suivant

          enddo

10        continue

        enddo

      Else

! .....On descend vers les fils

        if( associated(maille%fils(1)%ptr) ) then
          do l = 1, 2**ndim
            call feuille_fictive(maille%fils(l)%ptr)
          enddo
        endif

      Endif

      end subroutine feuille_fictive
