      subroutine support_racine(root)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      TYPE(structure_fils), &
        DIMENSION(nbre_racine(1), nbre_racine(2), nbre_racine(3)) :: root

      INTEGER :: ir, jr, kr
      INTEGER :: inc, k, l, m, n

      TYPE(structure_maille), POINTER :: courant

! =================== Allocation du support de la racine =============

      do kr = 1, nbre_racine(3)
        do jr = 1, nbre_racine(2)
          do ir = 1, nbre_racine(1)

            courant => root(ir,jr,kr)%ptr

            if( ir+1 <= nbre_racine(1) ) &
              courant%support(1)%suivant => root(ir+1,jr,kr)%ptr

            if( ir-1 >= 1 ) &
              courant%support(1)%precedent => root(ir-1,jr,kr)%ptr

          enddo

! .....Condition de periodicite
  
          if( caltype_deb(1) == 2 .or. caltype_fin(1) == 2 ) then

            root(1,jr,kr)%ptr%support(1)%precedent => &
                                            root(nbre_racine(1),jr,kr)%ptr

            root(nbre_racine(1),jr,kr)%ptr%support(1)%suivant => &
                                            root(1,jr,kr)%ptr

          endif

        enddo
      enddo

      if( ndim >= 2 ) then
       
         do kr = 1, nbre_racine(3)
          do ir = 1, nbre_racine(1)
            do jr = 1, nbre_racine(2)

              courant => root(ir,jr,kr)%ptr

              if( jr+1 <= nbre_racine(2) ) &
                courant%support(2)%suivant => root(ir,jr+1,kr)%ptr

              if( jr-1 >= 1 ) &
                courant%support(2)%precedent => root(ir,jr-1,kr)%ptr

            enddo

! .....Condition de periodicite
  
            if( caltype_deb(2) == 2 .or. caltype_fin(2) == 2 ) then

              root(ir,1,kr)%ptr%support(2)%precedent => &
                                 root(ir,nbre_racine(2),kr)%ptr

              root(ir,nbre_racine(2),kr)%ptr%support(2)%suivant => &
                                 root(ir,1,kr)%ptr

            endif

          enddo
        enddo

      endif

      if( ndim == 3 ) then
       
        do jr = 1, nbre_racine(2)
          do ir = 1, nbre_racine(1)
            do kr = 1, nbre_racine(3)

              courant => root(ir,jr,kr)%ptr

              if( kr+1 <= nbre_racine(3) ) &
                courant%support(3)%suivant => root(ir,jr,kr+1)%ptr

              if( kr-1 >= 1 ) &
                courant%support(3)%precedent => root(ir,jr,kr-1)%ptr

            enddo

! .....Condition de periodicite
  
            if( caltype_deb(3) == 2 .or. caltype_fin(3) == 2 ) then

              root(ir,jr,1)%ptr%support(3)%precedent => &
                                    root(ir,jr,nbre_racine(3))%ptr

              root(ir,jr,nbre_racine(3))%ptr%support(3)%suivant => &
                                    root(ir,jr,1)%ptr

            endif

          enddo
        enddo

      endif

      end subroutine support_racine
