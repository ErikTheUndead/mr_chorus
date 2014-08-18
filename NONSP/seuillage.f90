      recursive subroutine seuillage(maille)

      use mod_common
      use mod_structure
      use mod_recherche

      IMPLICIT NONE

      TYPE(structure_maille), pointer :: maille

      INTEGER :: i, j, k, l
      INTEGER :: indice, increment
      INTEGER :: niv_sup, sj, sk
      INTEGER, DIMENSION(ndim) :: ind_loc, ind_sup, ind_cor
      INTEGER, DIMENSION(3) :: ind

      TYPE(structure_maille), pointer :: courant

      REAL(KIND = 8) :: precis

! ============= seuillage et marquage des branches a elaguer =========

! *****C'est une feuille

      if( .not.associated(maille%fils(1)%ptr) ) then

! .....des fils ne doivent pas etre crees

        maille%arbre = .false.

! .....Mais si le pere d'une feuille reelle a un detail important 
!           ==> creation potentielle

        if( maille%feuille ) then        

          if((abs(maille%pere%detail / det_max) >= &
              maille%pere%eps) )  maille%arbre = .true.

        endif
        
      else

! *****On passe aux fils

        do l = 1, 2**ndim
          call seuillage(maille%fils(l)%ptr)
        enddo

! *****test sur le pere

! .....le detail est negligeable ==> point potentiellement plus dans arbre

        if((abs(maille%detail / det_max) < maille%eps) &
          .and. (maille%maille_niv > niv_min)) then

          maille%arbre = .false.

! .....sinon, le point est potentiellement dans arbre
        else
          maille%arbre = .true.

! .....niveau et indice de la maille courante
          niv_sup = maille%maille_niv
          ind_sup(:) = maille%maille_i(:)

! .....ainsi que ces proches voisins pour reconstruction polynomiale
!      si polynome de degre >= 1
          if(s > 0 ) then

            if(ndim == 3) then
              sj = s
              sk = s
            elseif(ndim == 2) then
              sj = s
              sk = 0
            elseif(ndim == 1) then
              sj = 0
              sk = 0
            endif

! .....Recherche de tous les points du support
            do k = -sk, sk
              ind(3) =  k
              do j = -sj, sj
                ind(2) = j
                do i = -s, s
                  ind(1) = i

                  do l = 1,ndim
                    ind_loc(l) = ind_sup(l) + ind(l)
                  enddo

! .....conditions aux limites
                  call cal_indice(niv_sup, ind_loc, ind_cor)

! .....recherche du point dans l'arbre 
                  indice = 1
                  increment = 1
                  do l = 1, ndim
                    indice = indice + (ind_cor(l) - 1) * increment
                    increment = increment * 2**niv_sup * nbre_racine(l)
                  enddo
                  indice = indice + (2**(ndim*niv_sup) - 1)/(2**ndim-1) * &
                            nbre_racine(1)*nbre_racine(2)*nbre_racine(3)

                  if( .not.associated(hash_table(indice)%ptr) ) then

                    print*, ' !! Seuillage !! : ', &
                            ' Pb de recherche support maille - ', &
                            'Maille niv = ',niv_sup,' Maille ijk = ',ind_sup, &
                            ' i, j, k loc = ',ind_loc,' i, j, k cor = ',ind_cor
                    stop
                  endif

                  courant => hash_table(indice)%ptr          

                  courant%arbre = .true.

                enddo
              enddo
            enddo

          endif

! +++++test sur croissance des singularites

! .....norme L1
          precis = 2.**(2*s-1) * maille%eps
! .....norme L2
!!          precis = 2.**((2*s-1)/2.) * maille%eps

          if(abs(maille%detail / det_max) >= precis ) then

            maille%arbre = .true.

            do l = 1, 2**ndim
              maille%fils(l)%ptr%arbre = .true.
            enddo

          endif

        endif

      endif

      end subroutine seuillage
