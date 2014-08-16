      recursive subroutine chainage(maille, Profondeur, Hauteur, &
                                      Position, Pere_maille)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      INTEGER :: i, k, l, m, n

      INTEGER, intent(in) :: Profondeur, Hauteur
      INTEGER, DIMENSION(ndim), intent(in) :: Position
      INTEGER, DIMENSION(ndim) :: Pos
      TYPE(structure_maille), POINTER :: maille
      TYPE(structure_maille), POINTER :: Pere_maille

      if (.not.associated(maille)) then
  
! .....On alloue la maille

        allocate(maille)

! +++++Valorisation de la position
 
        maille%maille_niv=Hauteur  
        maille%maille_i(:)=Position(:)

! +++++Valorisation de seuil suivant le niveau

! .....norme L1  
        maille%eps=epsilon * ( 2.**(ndim*(maille%maille_niv - niv_max)) )
! .....norme L2
!!        maille%eps=epsilon * ( 2.**(ndim*(maille%maille_niv - niv_max)/2.) )

! +++++Positionnement des variables logiques

! .....Les points appartenant a l'arbre seront initialises plus tard
        maille%arbre = .false.

! .....Les feuilles sont initialisees reelles (non fictives)
        maille%feuille = .true.

! +++++Valorisation du chainage vers le pere
 
        maille%pere => Pere_maille

! +++++Initialisation  des liens

        do l = 1, 2**ndim
          nullify(maille%fils(l)%ptr)
        enddo

        do l = 1,ndim
          nullify(maille%support(l)%precedent)
          nullify(maille%support(l)%suivant)
        enddo
  
      endif
    
      if( Profondeur > Hauteur) then

        do l = 1, 2**ndim

          k = 0
          n = 0
          do m = ndim, 1, -1
            k = k+1
            i = ((l-1) - n)/2**(ndim-k)  
            n = n + i*2**(ndim-k)

            Pos(m) = 2*Position(m)+i-1
          enddo

          call chainage(maille%fils(l)%ptr, Profondeur, Hauteur+1, &
                        Pos, maille)

        enddo
  
      endif

      end subroutine chainage
