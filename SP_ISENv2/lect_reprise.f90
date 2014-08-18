      subroutine lect_reprise(root)

      use mod_common
      use mod_structure
      use mod_arbre

      IMPLICIT NONE

      TYPE(structure_fils), &
         DIMENSION(nbre_racine(1), nbre_racine(2), nbre_racine(3)) :: root
 
      LOGICAL :: flag, info

      INTEGER :: ir, jr, kr, l, m, nbf

      INTEGER :: niv_lec
      INTEGER, DIMENSION(ndim) :: ijk_lec

      REAL(kind=8), DIMENSION(ndim) :: x_lec
      REAL(kind=8), DIMENSION(nvar) :: u_lec

      TYPE(structure_maille), POINTER :: root_loc
      TYPE(structure_maille), POINTER :: maille

! =================== Lecture des feuilles de l'arbre ==================

! +++++Nbre d'iterations, temps de simulation et nombre de feuilles

      read(20) nbre_iter, temps, compteur_feuille

      Print*
      Print*, ' Reprise d''une solution : Iterations = ', nbre_iter
      Print*, '                               Temps = ', temps
      Print*, '                    Nbre de feuilles = ', compteur_feuille
      Print*

! *****Pour toutes les feuilles de l'arbre

      do nbf = 1, compteur_feuille

        read(20, end = 10) ( ijk_lec(l), l=1,ndim ), &
                  niv_lec, &
                  ( x_lec(l), l=1,ndim ), &
                  ( u_lec(m), m=1,nvar )

! +++++Recherche de la maille pour validation des champs de la structure 

! .....recherche de la racine appropriee
        do kr = 1, nbre_racine(3)
          do jr = 1, nbre_racine(2)
            do ir = 1, nbre_racine(1)

              flag = .true.

              do l = 1, ndim
                if( ijk_lec(l) > &
          root(ir,jr,kr)%ptr%maille_i(l) * 2**niv_lec )  flag = .false.
              enddo
 
              if( flag ) goto 20

            enddo
          enddo
        enddo

        if( .not.flag ) then
          Print*, ' Lect_reprise : Probleme de localisation racine !! '
          Stop
        endif

20      continue

! .....recherche du point dans la racine appropriee
        root_loc => root(ir,jr,kr)%ptr

        call recherche(root_loc, ijk_lec, niv_lec, maille, info)

! +++++le point est trouve dans l'arbre
        if( info ) then

! .....valorisation des champs

          maille%maille_i(:) = ijk_lec(:)
          maille%maille_niv = niv_lec

          maille%x(:) = x_lec(:)
          maille%u(:) = u_lec(:)

! .....elagage du sous arbre

          call elagage(maille)

! +++++le point n'est pas dans l'arbre !!
        else

          Print*, ' Lect_reprise : Probleme !! '
          Print*, '                => le point n''est pas dans l''arbre ' 
          Print*, ' Point recherche : niv = ', niv_lec, &
                  '                   i,j,k = ', ijk_lec
         
          Stop

        endif

      enddo

! +++++Fin de lecture
      return

! +++++Probleme de lecture
10    continue

      Print*, ' Probleme de lecture du fichier reprise !! '
      Print*, '    - Nbre de points lus = ', nbf-1
      Stop

      end subroutine lect_reprise
