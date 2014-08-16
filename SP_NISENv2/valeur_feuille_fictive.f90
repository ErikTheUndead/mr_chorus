      recursive subroutine valeur_feuille_fictive(maille)

      use mod_common
      use mod_structure
      use mod_recherche

      implicit none

      TYPE(structure_maille), pointer :: maille

      LOGICAL :: flag, reelle

      INTEGER :: l, m, s_deb, st

      TYPE(structure_maille), pointer :: courant

! +++++On valorise les feuilles fictives a partir des peres

      If( associated( maille%fils(1)%ptr) ) then

        if( .not.maille%fils(1)%ptr%feuille ) then

          courant => maille

! *****Valorisation des fils fictifs

          call decoding_moy(courant)

        endif

! .....On descend vers les fils

          do l = 1, 2**ndim
            call valeur_feuille_fictive(maille%fils(l)%ptr)
          enddo
  
      Endif

      end subroutine valeur_feuille_fictive
