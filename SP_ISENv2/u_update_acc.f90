	recursive subroutine u_update_acc(maille)

	use mod_common
	use mod_structure

	implicit none

        INTEGER :: l, m

	TYPE(structure_maille), pointer :: maille

! +++++on parcours l'arbre a partir de la racine vers les feuilles

        If( associated(maille) ) then

! .....On initialise la solution aux temps intermediaires

          do m = 2,nvar
            maille%u(m) = maille%u(m)/maille%u(1)
          end do

! +++++On passe au fils pour la mise a jour

          if( associated(maille%fils(1)%ptr) ) then
            do l = 1, 2**ndim
              call u_update_acc(maille%fils(l)%ptr)
            enddo
          endif
        Endif

	end subroutine u_update_acc