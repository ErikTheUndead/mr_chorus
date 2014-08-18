      subroutine creation(maille, reelle)

      use mod_common
      use mod_structure

      IMPLICIT NONE	

      LOGICAL :: reelle

      TYPE(structure_maille), pointer :: maille

      INTEGER :: i, k, l, lm, m, n
      INTEGER :: indice, increment

      REAL(kind=8) :: delta

! +++++On alloue et on valorise les fils

      do l = 1, 2**ndim

! *****Allocation de maille_fils

        allocate(maille%fils(l)%ptr)

! *****valorisation des valeurs de la structure

! .....niveau
        maille%fils(l)%ptr%maille_niv= maille%maille_niv+1		

! .....indices
        k = 0
        n = 0
        do m = ndim, 1, -1
          k = k+1
          i = ((l-1) - n)/2**(ndim-k)
          n = n + i*2**(ndim-k)

          maille%fils(l)%ptr%maille_i(m) = 2*maille%maille_i(m)+i-1

! .....coordonnees
          delta = ( float(i)*.5 - 0.25 ) * maille%dx(m)
          maille%fils(l)%ptr%x(m)= maille%x(m) + delta

! .....pas d'espace
          maille%fils(l)%ptr%dx(m)= .5*maille%dx(m)

        enddo
        
! .....la maille a-t-elle des fils ?
        maille%fils(l)%ptr%arbre = .false.
        
! .....la feuille est-elle reelle ou fictive
        maille%fils(l)%ptr%feuille = reelle

! .....niveau d'erreur admissible
! .....norme L1
        maille%fils(l)%ptr%eps= epsilon * &
                               (2.**(ndim*(maille%maille_niv-niv_max)))
! .....norme L2
!!        maille%fils(l)%ptr%eps= epsilon * &
!!                               (2.**(ndim*(maille%maille_niv-niv_max)-1))

! .....details
        do lm = 1, 2**ndim
          do m = 1,nvar
            maille%fils(l)%ptr%det(m, lm) = 0.
          enddo
        enddo

        maille%fils(l)%ptr%detail = 0.

! .....initialisation des flux

        maille%fils(l)%ptr%fmp_dt(:) = 0.
        maille%fils(l)%ptr%fmp_gch(:) = 0.

        maille%fils(l)%ptr%fvisc_dt(:,:) = 0.
        maille%fils(l)%ptr%fvisc_gch(:,:) = 0.
                
! *****chainage	

        maille%fils(l)%ptr%pere => maille		

! *****Entree dans la table de hashage

! .....reconstruction de l'indice dans la table
        indice = 1 
        increment = 1
        do lm = 1, ndim
          indice = indice + (maille%fils(l)%ptr%maille_i(lm) - 1) * increment
 
          increment = increment * 2**maille%fils(l)%ptr%maille_niv * &
                                     nbre_racine(lm)
        enddo
        indice = indice + (2**(ndim*maille%fils(l)%ptr%maille_niv) - 1) / &
                 (2**ndim-1) * nbre_racine(1)*nbre_racine(2)*nbre_racine(3)

! .....table de hashage
        hash_table(indice)%ptr => maille%fils(l)%ptr

      enddo

      end subroutine creation
