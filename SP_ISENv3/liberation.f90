      subroutine liberation(maille)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      INTEGER :: l, m
      INTEGER :: indice, increment

      TYPE(structure_maille), pointer :: maille

      TYPE(structure_maille), pointer :: courant

! +++++On desalloue les fils et toute la structure associee

      do l = 1, 2**ndim

        if( associated(maille%fils(l)%ptr) ) then

          courant => maille%fils(l)%ptr

! .....table de hashage ==> null
          indice = 1 
          increment = 1
          do m = 1, ndim
            indice = indice + (courant%maille_i(m) - 1) * increment
            increment = increment * 2**courant%maille_niv * nbre_racine(m)
          enddo
          indice = indice + (2**(ndim*courant%maille_niv) - 1)/(2**ndim-1) * &
                            nbre_racine(1)*nbre_racine(2)*nbre_racine(3)

! .....table de hashage
          nullify( hash_table(indice)%ptr )

          deallocate(maille%fils(l)%ptr)

        endif

      enddo

      end subroutine liberation
