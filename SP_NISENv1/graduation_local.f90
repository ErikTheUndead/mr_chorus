      recursive subroutine graduation_local(maille, reelle)

      use mod_common
      use mod_structure
      use mod_recherche

      IMPLICIT NONE

      TYPE(structure_maille), pointer :: maille

      LOGICAL :: reelle

      INTEGER :: i, j, k, l
      INTEGER :: indice, increment
      INTEGER :: niv_sup, niv_pere, sj, sk
      INTEGER, DIMENSION(ndim) :: ind_loc, ind_sup, ind_crt
      INTEGER, DIMENSION(3) :: ind

      TYPE(structure_maille), pointer :: courant

! =========== recherche du support et creation si necessaire =========

      if(ndim == 3) then
        sj = sgrad
        sk = sgrad
      elseif(ndim == 2) then
        sj = sgrad
        sk = 0
      elseif(ndim == 1) then
        sj = 0
        sk = 0
      endif

! .....niveau et indice de la maille courante
      niv_sup = maille%maille_niv
      ind_sup(:) = maille%maille_i(:)

! .....Recherche de tous les points du support
      do k = -sk, sk
        ind(3) =  k
        do j = -sj, sj
          ind(2) = j
          do i = -sgrad, sgrad
            ind(1) = i

            do l = 1,ndim
              ind_loc(l) = ind_sup(l) + ind(l)
            enddo

! .....conditions aux limites
            call cal_indice(niv_sup, ind_loc, ind_crt)

! .....recherche du point dans l'arbre 
            indice = 1
            increment = 1
            do l = 1, ndim
              indice = indice + (ind_crt(l) - 1) * increment
              increment = increment * 2**niv_sup * nbre_racine(l)
            enddo
            indice = indice + (2**(ndim*niv_sup) - 1)/(2**ndim-1) * &
                              nbre_racine(1)*nbre_racine(2)*nbre_racine(3)

            if( associated(hash_table(indice)%ptr) ) then

              courant => hash_table(indice)%ptr          

              courant%pere%arbre = .true.

! .....si les fils sont dans l'arbre ==> le point doit etre dans l'arbre
              if(associated(courant%fils(1)%ptr)) then
                do l = 1, 2**ndim

                  if( (courant%fils(l)%ptr%arbre) .and. &
                      (.not.courant%arbre) ) courant%arbre = .true.

                enddo
              endif

! .....le point n'est pas dans l'arbre ==> il doit etre creer

            else

! .....indice et niveau du pere
              niv_pere = niv_sup - 1
              ind_loc(:) = (ind_crt(:) -1)/2 + 1

! .....conditions aux limites
              call cal_indice(niv_pere, ind_loc, ind_crt)

! .....recherche du point dans l'arbre 
              indice = 1
              increment = 1
              do l = 1, ndim
                indice = indice + (ind_crt(l) - 1) * increment
                increment = increment * 2**niv_pere * nbre_racine(l)
              enddo
              indice = indice + (2**(ndim*niv_pere) - 1)/(2**ndim-1) * &
                                nbre_racine(1)*nbre_racine(2)*nbre_racine(3)

              if( .not.associated(hash_table(indice)%ptr) ) then

                Print*, ' Le Pere non plus n''est pas dans l''arbre ',&
                        ' ==> Creation !! '
                Print*, ' Pere : niv = ', niv_pere, ' ijk = ', ind_crt
                Stop
              endif

! .....Creation des fils
              courant => hash_table(indice)%ptr

              call creation(courant, reelle)

! .....son pere doit etre dans l'arbre
              courant%arbre = .true.

              compteur_create = compteur_create + 1

! .....Graduation
              call graduation_local(courant, reelle)

! .....Valorisation des Fils
              call decoding_moy(courant)

            endif

          enddo
        enddo
      enddo
!
      end subroutine graduation_local
