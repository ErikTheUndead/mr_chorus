      recursive subroutine graduation(maille)

      use mod_common
      use mod_structure
      use mod_fonction
      use mod_recherche

      IMPLICIT NONE

      TYPE(structure_maille), POINTER :: maille

      LOGICAL :: reelle

      INTEGER :: l

! ======== Graduation en partant de la racine vers les feuilles ======

      if( .not.associated(maille%fils(1)%ptr) ) then

! ++++++Feuille avec arbre = true ==> creation des fils si niveau < niv_max

        if((maille%arbre) .and. &
           (maille%maille_niv < niv_max) ) then

          reelle = .true.
          call creation(maille, reelle)
          compteur_create = compteur_create + 1

! .....Graduation
          call graduation_local(maille, reelle)

! .....Valorisation des Fils
          call decoding_moy(maille)

        endif

! +++++Sinon, on traite le pere

      else

! +++++Verification de la gradualite et creation si necessaire

        if((maille%arbre) .and. (maille%maille_niv >= niv_min)) then
          reelle = .true.
          call graduation_local(maille, reelle)
        endif

! +++++On traite les fils ensuite

        do l = 1, 2**ndim
          call graduation(maille%fils(l)%ptr)
        enddo

      endif

      end subroutine graduation
