      recursive subroutine init_flags(maille)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      INTEGER :: l, indice, increment
      TYPE(structure_maille), pointer :: maille

! +++++Initialisation de la table de hashage

      If ( associated(maille) ) then

! .....reconstruction de l'indice dans la table
        indice = 1 
        increment = 1
        do l = 1, ndim
          indice = indice + (maille%maille_i(l) - 1) * increment
 
          increment = increment * 2**maille%maille_niv * nbre_racine(l)
        enddo
        indice = indice + (2**(ndim*maille%maille_niv) - 1) / &
                 (2**ndim-1) * nbre_racine(1)*nbre_racine(2)*nbre_racine(3)

! .....table de hashage
        hash_table(indice)%ptr => maille

! +++++Initialisation des mailles dans l'arbre

        maille%arbre = .false.

! .....On descend vers les fils
        if(associated(maille%fils(1)%ptr)) then

          do l = 1, 2**ndim
            call init_flags(maille%fils(l)%ptr)
          enddo

        endif

      Endif

      end subroutine init_flags
