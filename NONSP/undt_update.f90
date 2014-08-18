        recursive subroutine undt_update(maille)

	use mod_common
	use mod_structure

	implicit none
 
        INTEGER :: l, m

	TYPE(structure_maille), pointer :: maille

! +++++on parcours l'arbre a partir de la racine vers les feuilles

        If( associated(maille) ) then

! .....On initialise la solution au temps n.Dt

          do m = 1,nvar
            maille%u_ndt(m) = maille%u(m)
          enddo

! +++++On passe au fils pour la mise a jour

          if( associated(maille%fils(1)%ptr) ) then
            do l = 1, 2**ndim
              call undt_update(maille%fils(l)%ptr)
            enddo
          endif
        Endif

	end subroutine undt_update
