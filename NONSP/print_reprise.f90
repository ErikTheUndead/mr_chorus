      subroutine print_reprise

      use mod_common
      use mod_structure

      IMPLICIT NONE
   
      INTEGER :: l, leaf, m

      TYPE(structure_maille), pointer :: maille

! +++++On parcourt les feuilles reelles (non-fictives)

      Do leaf = 1, compteur_feuille

        maille => feuille_reelle(leaf)%ptr

! +++++On ecrit les feuilles reelles (non-fictives)

        write(20) ( maille%maille_i(l), l=1,ndim ), &
                    maille%maille_niv, &
                  ( maille%x(l), l=1,ndim ), &
                  ( maille%u(m), m=1,nvar )

      Enddo
    
      end subroutine print_reprise
