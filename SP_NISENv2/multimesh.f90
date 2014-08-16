! .....subroutine qui genere un maillage regulier avec
! .....des niveaux de grilles niv allant de 0 (gros) a niv_max (fin)
! .....definit les limites physiques du domaine de calcul

      recursive subroutine multimesh(maille, x_debut, x_fin)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      INTEGER :: i, k, l, m, n

      TYPE(structure_maille), pointer :: maille
      REAL(kind=8), DIMENSION(ndim) :: x, x_debut, x_fin, dx
      REAL(kind=8), DIMENSION(ndim) :: x_beg, x_end

!
! +++++distribution des points du maillage et limites du domaine
!	
      if (associated(maille)) then

        dx(:)=x_fin(:)-x_debut(:)
        x(:)=x_debut(:)+dx(:)*.5

        do l = 1,2**ndim

          k = 0
          n = 0
          do m = ndim, 1, -1
            k = k+1
            i = ((l-1) - n)/2**(ndim-k)
            n = n + i*2**(ndim-k)

            x_beg(m) = x_debut(m)
            x_end(m) = x(m)
            if(i /= 0) then
              x_beg(m) = x(m)
              x_end(m) = x_fin(m)
            endif
          enddo

          call multimesh(maille%fils(l)%ptr, x_beg, x_end)

        enddo

        maille%x(:)=x(:)
        maille%dx(:)=dx(:)

      else 
        return
      endif  
   
      end subroutine multimesh
