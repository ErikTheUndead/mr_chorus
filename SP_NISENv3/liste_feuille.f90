      recursive subroutine liste_feuille(maille, compt)

      use mod_common
      use mod_structure

      INTEGER :: l
      INTEGER :: compt
      TYPE(structure_maille), POINTER :: maille

! .....parcours de l'arbre en postfix : fils ensuite pere

      If( .not.associated(maille%fils(1)%ptr) ) then

! .....C'est une feuille => mise dans la liste des vraies feuilles

        maille%feuille = .true.

        compt=compt+1

        feuille_reelle(compt)%ptr => maille

      Else

! .....On remonte vers les peres

        do l = 1, 2**ndim
          call liste_feuille(maille%fils(l)%ptr, compt)
        enddo

! .....on marque le noeud comme reel
        maille%feuille = .true.

      Endif

      end subroutine liste_feuille
