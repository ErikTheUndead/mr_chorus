!	recherche par indices (i,niv) <=> (ind, niv)

      recursive subroutine recherche(courant, ind, niv_loc, val_out, info)

      use mod_common
      use mod_structure

      implicit none

      type(structure_maille), pointer :: courant
      type(structure_maille), pointer :: val_out

      logical :: info, flag
      integer :: l, m
      integer :: niv_loc
      integer, dimension(ndim) :: ind, ind_loc

      info = .false.

      if(.not.associated(courant))then

        info = .false.
        nullify(val_out)
        return
  
      else  

        flag = .false.
        if (courant%maille_niv==niv_loc )  then
          do l = 1, ndim
            if(courant%maille_i(l) == ind(l)) then
              flag = .true.
            else
              flag = .false.
            endif

            if( .not.flag) exit
          enddo
        endif
  
        if ( flag )  then
    
          val_out => courant
          info = .true.
          return
    
        else
!
          l = 1
          do m = ndim, 1, -1
            ind_loc(m)=(2*courant%maille_i(m)-1)* &
                        2**(niv_loc-courant%maille_niv-1)

            if(ind(m) > ind_loc(m)) l = l + 2**(m-1) 

          enddo

          call recherche(courant%fils(l)%ptr, ind, niv_loc, val_out, info)
          if(info) return 
        endif  
      endif
 
      end subroutine recherche
