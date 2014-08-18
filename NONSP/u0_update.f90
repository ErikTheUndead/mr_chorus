	recursive subroutine u0_update(maille)

	use mod_common
	use mod_structure

	implicit none
 
        INTEGER :: l, m

	TYPE(structure_maille), pointer :: maille

! +++++on parcours l'arbre a partir de la racine vers les feuilles

        If( associated(maille) ) then

! .....On initialise la solution aux temps intermediaires
          do m = 1,nvar
            maille%u0(m) = maille%u(m)
          enddo

! +++++On passe au fils pour la mise a jour

          if( associated(maille%fils(1)%ptr) ) then
            do l = 1, 2**ndim
              call u0_update(maille%fils(l)%ptr)
            enddo
          endif

        Endif

	end subroutine u0_update
