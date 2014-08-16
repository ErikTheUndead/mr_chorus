      recursive subroutine elagage(maille)

      use mod_common
      use mod_structure
      use mod_fonction

      IMPLICIT NONE

      INTEGER :: l, m

      TYPE(structure_maille), pointer :: maille

! =============== Parcours de l'arbre a partir des feuilles ==============

      if( associated(maille%fils(1)%ptr) ) then

! +++++On passe aux fils

        do l = 1, 2**ndim
          call elagage(maille%fils(l)%ptr)
        enddo

! +++++On supprime les fils qui ne sont plus utiles

        if(.not.maille%arbre) then

          call liberation(maille)

        else
          maille%arbre = .true.
        endif

      endif

      end subroutine elagage
