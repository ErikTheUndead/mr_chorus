      subroutine mise_a_jour

      use mod_common
      use mod_structure

      implicit none

      INTEGER :: l, leaf, m

      TYPE(structure_maille), pointer :: maille


! +++++On parcourt les feuilles reelles (non-fictives)

      Do leaf = 1, compteur_feuille

        maille => feuille_reelle(leaf)%ptr

! .....On initialise la solution au temps n.Dt

        maille%u(:) = maille%unp1(:)

      Enddo

      end subroutine mise_a_jour
