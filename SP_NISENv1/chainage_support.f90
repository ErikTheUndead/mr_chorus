      recursive subroutine chainage_support(Pere_maille)

      use mod_common
      use mod_structure
	
      IMPLICIT NONE

      INTEGER :: inc, k, l, m, n
      INTEGER, DIMENSION(ndim) :: ijk

      TYPE(structure_maille), POINTER :: Pere_maille
      TYPE(structure_maille), POINTER :: maille

      TYPE(structure_maille), POINTER :: courant

!
! ========== Parcours de l'arbre a partir de la racine =============
!

      if( associated(Pere_maille%fils(1)%ptr) ) then

! +++++Allocation du support des fils

        do l = 1, 2**ndim

! .....Pour chaque fils
          maille => Pere_maille%fils(l)%ptr

! .....position de la maille
          k = 0
          n = 0
          do m = ndim, 1, -1
            k = k+1
            ijk(m) = ((l-1) - n)/2**(ndim-k)
            n = n + ijk(m)*2**(ndim-k)
          enddo

! .....Allocation et Valorisation du support pour la maille courante

          do m = 1, ndim
            inc = 2**(m-1)

! *****Pour le fils gauche

            if( ijk(m) == 0 ) then

! .....le suivant est son frere droit
              maille%support(m)%suivant => Pere_maille%fils(l+inc)%ptr

! .....le precedent est son cousin droit
              if( associated(Pere_maille%support(m)%precedent) ) then
                courant => Pere_maille%support(m)%precedent

                if( associated(courant%fils(l+inc)%ptr) ) then

                  maille%support(m)%precedent => courant%fils(l+inc)%ptr

                else
                  
                  nullify(maille%support(m)%precedent)

                endif

              else
                  
                nullify(maille%support(m)%precedent)

              endif

! *****Pour le fils droit

            else

! .....le precedent est son frere gauche
              maille%support(m)%precedent => Pere_maille%fils(l-inc)%ptr

! .....le suivant est son cousin gauche
              if(associated(Pere_maille%support(m)%suivant) ) then
                courant => Pere_maille%support(m)%suivant

                if( associated(courant%fils(l-inc)%ptr) ) then

                  maille%support(m)%suivant => courant%fils(l-inc)%ptr

                else

                  nullify(maille%support(m)%suivant)

                endif

              else

                nullify(maille%support(m)%suivant)

              endif

            endif

          enddo

        enddo

        do l = 1, 2**ndim
          call chainage_support(Pere_maille%fils(l)%ptr)
        enddo
 
      endif

      end subroutine chainage_support
