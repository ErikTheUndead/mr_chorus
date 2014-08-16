	recursive subroutine encoding_moy(maille)
	
	use mod_common
	use mod_structure
	
	IMPLICIT NONE

        INTEGER :: l, m

	TYPE(structure_maille), pointer :: maille

! .....parcours de l'arbre en postfix : fils ensuite pere

	if( .not.associated(maille%fils(1)%ptr) ) then
	
! feuilles donc on ne fait rien	
	
	else

! .....Dans le cas ou la feuille est reelle (pas fictive)
 
          if( maille%fils(1)%ptr%feuille ) then

! .....valorisation des peres
	
            do m = 1,nvar
	      maille%u(m) = 0.
            enddo

            do l = 1, 2**ndim
	  
	      call encoding_moy(maille%fils(l)%ptr)		
!
              do m = 1,nvar
	        maille%u(m) = maille%u(m) + &
                              (maille%fils(l)%ptr%u(m)) / float(2**ndim)
              enddo
		
            enddo

	  endif	

	endif

	end subroutine encoding_moy
