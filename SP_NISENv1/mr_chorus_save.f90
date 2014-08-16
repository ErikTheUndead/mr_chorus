      Program MR_CHORUS

      use mod_common
      use mod_structure
      use mod_initial
      use mod_arbre
      use mod_recherche
      use mod_print
      use mod_fonction
      use mod_euler
      use mod_save_vtk

      IMPLICIT NONE

      LOGICAL :: cok

      CHARACTER(len=80) :: nom_fich

      INTEGER :: i, j, k, i_sp
      INTEGER :: dir, l, m, n_st, n_ite, strang_pas

      INTEGER, DIMENSION(3,6) :: strang
      INTEGER, DIMENSION(3) :: ind_pos
      INTEGER, DIMENSION(ndim) :: position

      REAL(kind=8) :: err_Etot
      REAL(kind=8), DIMENSION(ndim) :: x_deb, x_fin
      REAL(kind=8) :: min_ang
      REAL(kind=8) :: TempsDeb,TempsFin


      TYPE(structure_maille), Pointer :: root, courant, work
      TYPE(structure_fils), Pointer, DIMENSION(:,:,:) :: racine

! ====================== Initialisation ============================

! .....definition du nombre Pi

      pi=acos(-1.)

! .....initialisation de la clef de fin de calcul
      cok = .false.

! +++++definition de la tensorisation pour symetrie

      strang(:,:) = 0.

      if(ndim == 1) then

        strang_pas = 1

        strang(1,1) = 1

      elseif(ndim == 2) then

        strang_pas = 2

        strang(1,1) = 1
        strang(2,1) = 2

        strang(1,2) = 2
        strang(2,2) = 1

      elseif(ndim == 3) then

        strang_pas = 6

        strang(1,1) = 1
        strang(2,1) = 2
        strang(3,1) = 3

        strang(1,2) = 1
        strang(2,2) = 3
        strang(3,2) = 2

        strang(1,3) = 2
        strang(2,3) = 3
        strang(3,3) = 1

        strang(1,4) = 2
        strang(2,4) = 1
        strang(3,4) = 3

        strang(1,5) = 3
        strang(2,5) = 2
        strang(3,5) = 1

        strang(1,6) = 3
        strang(2,6) = 1
        strang(3,6) = 2

      endif

! *****Ouverture fichier pour ecriture des residus

      open(10, file='res_mrchorus.dat',form='formatted')

! ======================= lecture des donnees =========================

      call lect_don


      write(nom_fich,995) niv_max,s,epsilon

      open(15, file=nom_fich,form='formatted')

      write(15,*) 'TITLE=Choc_Spot_Ma=1.1588'
      write(15,*) 'VARIABLES="temps","Int_ w","Int_|w|","Int_GPGrho"', &
                  '"Int_|GPGrho|"'

! ========================== Initialisations ==========================

! .....definition du symbol de Kroenecker
      delta_ij(:,:) = 0.
      do m = 1, ndim
        delta_ij(m,m) = 1.
      enddo

! .....etendu du support pour la graduation
      sgrad = s + (nordre+1)/2
!!      sgrad = s + (nordre)/2 + 1

      Print*, ' '
      Print*, ' Sgrad = ', sgrad

! .....pas d'espace sur le niveau le plus fin
      do l = 1,ndim
        dxmin(l) = ( lim_fin(l) - lim_deb(l) ) / &
                   ( 2.**niv_max * nbre_racine(l) )
      enddo

! +++++Initialisation des coefficients d'interpolation selon l'ordre

! .....Allocation tableaux des coefficients d'interpolation d'ordre 2*s
      ALLOCATE(coef(s))

! .....ordre s = 1
      if(s == 1) then

        coef(1)=-0.125

! .....ordre s = 2
      elseif(s == 2) then

        coef(1)=-22./128.
        coef(2)=3./128.

! .....ordre s = 3
      elseif(s == 3) then

        coef(1)=-201/1024.
        coef(2)=11./256.
        coef(3)=-5./1024.

      endif

! +++++Initialisation des racines de l'arbre

! .....Allocation des racines
      ALLOCATE( racine(nbre_racine(1),nbre_racine(2), nbre_racine(3)) )

! .....Declaration memoire de la table de correspondance

      ALLOCATE( hash_table(Nbre_pt_arbre) )
      ALLOCATE( feuille_reelle(nbre_maille) )

! .....definition des indices maxima et du nombre de variables
      i_max = nbre_racine(1)*2**niv_max
      if( ndim == 1 ) then
        j_max = 1
        k_max = 1
        ijk_max = (i_max+1)
        nprint_var = nvar+3
      elseif( ndim == 2 ) then
        j_max = nbre_racine(2)*2**niv_max
        k_max = 1
        ijk_max = (i_max+1)*(j_max+1)
        nprint_var = nvar+2
      else
        j_max = nbre_racine(2)*2**niv_max
        k_max = nbre_racine(3)*2**niv_max
        ijk_max = (i_max+1)*(j_max+1)*(k_max+1)
        nprint_var = nvar+ndim+2
      endif

! ============== Definition de l'arbre et de la solution ==============

! +++++Creation et chainage de l'arbre

      do k = 1, nbre_racine(3)
        ind_pos(3) = k
        do j = 1, nbre_racine(2)
          ind_pos(2) = j
          do i = 1, nbre_racine(1)
            ind_pos(1) = i

            do l = 1, ndim
              position(l) = ind_pos(l)
            enddo
            call chainage(racine(i,j,k)%ptr, niv_max, 0, position, courant)

          enddo
        enddo
      enddo

! +++++Creation et chainage du support

! .....creation du support pour la racine

      call support_racine(racine)

! .....support pour tous l'arbre

      do k = 1, nbre_racine(3)
        do j = 1, nbre_racine(2)
          do i = 1, nbre_racine(1)

            root => racine(i, j, k)%ptr
            call chainage_support(root)

           enddo
        enddo
      enddo

! ++++++Definition du maillage

      do k = 1, nbre_racine(3)
        ind_pos(3) = k
        do j = 1, nbre_racine(2)
          ind_pos(2) = j
          do i = 1, nbre_racine(1)
            ind_pos(1) = i

            do l = 1, ndim
              x_deb(l) = lim_deb(l) + (lim_fin(l) - lim_deb(l)) * &
                         (ind_pos(l) - 1) / nbre_racine(l)

              x_fin(l) = lim_deb(l) + (lim_fin(l) - lim_deb(l)) * &
                         ind_pos(l) / nbre_racine(l)
            enddo
            root => racine(i, j, k)%ptr
            call multimesh(root, x_deb, x_fin)

          enddo
        enddo
      enddo

! +++++Initialisation de la solution sur le niveau le plus fin

! .....Initialisation du temps de simulation
        temps=0.
        nbre_iter = 0

! .....initialisation de la solution par valeurs moyennes u(i,niv_max)
!      sur la grille la plus fine

        do k = 1, nbre_racine(3)
          do j = 1, nbre_racine(2)
            do i = 1, nbre_racine(1)

              root => racine(i, j, k)%ptr
              call init_sol(root)

            enddo
          enddo
        enddo

! +++++Si reprise, lecture d'une solution sur un arbre existant

      if( cle_reprise ) then

        open(20, file='save_mrchorus.dat', form='unformatted')

        call lect_reprise(racine)

        close(20)

      endif

! .............................................................
! .......................Boucle en temps.......................
! .............................................................

	call cpu_time(TempsDeb)

      Do nt = 1, npdt

! .....intialisation de arbre = .false.

        do k = 1, nbre_racine(3)
          do j = 1, nbre_racine(2)
            do i = 1, nbre_racine(1)

              root => racine(i, j, k)%ptr

              call init_flags(root)

            enddo
          enddo
        enddo

! +++++Codage par valeurs moyennes
!      u connue par valeurs moyennes u(i,niv_max) sur grille la plus fine

        do k = 1, nbre_racine(3)
          do j = 1, nbre_racine(2)
            do i = 1, nbre_racine(1)

              root => racine(i, j, k)%ptr

              call encoding_moy(root)

            enddo
          enddo
        enddo

! +++++Calcul des valeurs de details sur tout l'arbre

        det_max = 0.

        do k = 1, nbre_racine(3)
          do j = 1, nbre_racine(2)
            do i = 1, nbre_racine(1)

              root => racine(i, j, k)%ptr

              call encoding_det(root)

            enddo
          enddo
        enddo

! ============== Generation de la grille hybride initiale ============

! +++++Construction de l'arbre hybride

        do k = 1, nbre_racine(3)
          do j = 1, nbre_racine(2)
            do i = 1, nbre_racine(1)

              root => racine(i, j, k)%ptr

              call seuillage(root)

            enddo
          enddo
        enddo

! +++++Construction de l'arbre hybride

        do
          compteur_create = 0

! .....verification de la graduation et creation si necessaire

          do k = 1, nbre_racine(3)
            do j = 1, nbre_racine(2)
              do i = 1, nbre_racine(1)

                root => racine(i, j, k)%ptr

                call graduation(root)

              enddo
            enddo
          enddo

! .....valorisation ==> decodage par valeur moyenne

          do k = 1, nbre_racine(3)
            do j = 1, nbre_racine(2)
              do i = 1, nbre_racine(1)

                root => racine(i, j, k)%ptr

                call decoding_moy(root)

              enddo
            enddo
          enddo


          if( compteur_create == 0 ) exit
        enddo

! +++++Elagage de l'arbre hybride

        if( .not.cle_reprise .or. nt /= 1 ) then

          do k = 1, nbre_racine(3)
            do j = 1, nbre_racine(2)
              do i = 1, nbre_racine(1)

                racine(i, j, k)%ptr%arbre = .true.

                root => racine(i, j, k)%ptr

                call elagage(root)

              enddo
            enddo
          enddo

        endif

! +++++Etablissament du support sur le nouvel arbre

        do k = 1, nbre_racine(3)
          do j = 1, nbre_racine(2)
            do i = 1, nbre_racine(1)

              root => racine(i, j, k)%ptr
              call chainage_support(root)

            enddo
          enddo
        enddo

! +++++Etablissement de la liste des feuilles

        compteur_feuille = 0

        do k = 1, nbre_racine(3)
          do j = 1, nbre_racine(2)
            do i = 1, nbre_racine(1)

              root => racine(i, j, k)%ptr
              call liste_feuille(root, compteur_feuille)

            enddo
          enddo
        enddo

! +++++Etablissament des feuilles fictives : molecule Euler

        compteur_fictive = 0

        do k = 1, nbre_racine(3)
          do j = 1, nbre_racine(2)
            do i = 1, nbre_racine(1)

              root => racine(i, j, k)%ptr

              call feuille_fictive(root)

            enddo
          enddo
        enddo

! +++++Etablissament du support sur le nouvel arbre
!
        do k = 1, nbre_racine(3)
          do j = 1, nbre_racine(2)
            do i = 1, nbre_racine(1)

              root => racine(i, j, k)%ptr
              call chainage_support(root)

            enddo
          enddo
        enddo

! +++++Etablissament des feuilles fictives : molecule visqueuse

        if( ndim > 1 .and. Reynolds /= 0. ) then

          do k = 1, nbre_racine(3)
            do j = 1, nbre_racine(2)
              do i = 1, nbre_racine(1)

                root => racine(i, j, k)%ptr
  !
                call feuille_fictive_visc(root)

              enddo
            enddo
          enddo

! +++++Etablissament du support sur le nouvel arbre

          do k = 1, nbre_racine(3)
            do j = 1, nbre_racine(2)
              do i = 1, nbre_racine(1)

                root => racine(i, j, k)%ptr
                call chainage_support(root)

              enddo
            enddo
          enddo

        endif

! +++++Ecriture de la solution initiale

        if( .not.cle_reprise .and. nt == 1 ) then



! .....Allocation des tableaux de travail pour l'ecriture
          ALLOCATE( flag_coins(ijk_max), x_coins(ijk_max,ndim) )

! .....allocation des tableaux coordonnees et champs
          ALLOCATE( ind_coins(compteur_feuille, 2**ndim) )
          ALLOCATE( var_print(compteur_feuille, nprint_var) )

          ALLOCATE( coords(compteur_feuille, ndim) )

! .....Fichier coordonnees des feuilles de l'arbre

          nbre_coins = 0

          flag_coins(:) = -1

          call print_leaves

          if( ndim == 1 ) then

! .....Fichier des champs sur grille complete
            write(nom_fich,999) int(100.*temps),niv_max,s,epsilon
            open(31,file=nom_fich,status='unknown')

! .....Ecriture au format colonne pour Gnuplot
            write(31,*) '# TITLE="1D_Dyadic_Tree"'
            write(31,*) '# VARIABLES="x","dx","niv","rho","u_x","P","T"'

            do l = 1, compteur_feuille
              write(31,900) x_coins(l,1),(var_print(l,m), m = 1, nprint_var)
            enddo

            close(31)


! .....Ecriture au format E-F pour Tecplot ou Paraview
          elseif( ndim == 2 ) then

!             write(31,*) 'TITLE="2D_Dyadic_Tree_OS7DIF2"'
!             write(31,*) 'VARIABLES="x","y","rho",', &
!                           '"u_x","u_y","P","T","w"'
!             write(31,*) 'ZONE T="Pure_QUAD"', &
!                     ', NODES = ', nbre_coins, &
!                     ', ELEMENTS = ', compteur_feuille, &
!                     ', DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL', &
!                     ', VARLOCATION=([1-2]=NODAL, [3-8]=CELLCENTERED)'
!
! ! .....ecriture de la liste des coins des cellules
!             write(31, 900)  ((x_coins(l, m), l = 1, nbre_coins), &
!                               m = 1, ndim )
!
! ! .....ecriture des variables au centre des cellules
!             write(31, 900) ((var_print(l, m), l = 1, compteur_feuille), &
!                              m = 1, nprint_var)
!
! ! .....ecriture de la liste de connectivite
!             do l = 1, compteur_feuille
!               write(31, 800)  (ind_coins(l, m), m = 1, 2**ndim)
!             enddo

            call save_VTK_2D

! .....Ecriture au format E-F pour Tecplot ou Paraview
          else

! .....Fichier des champs sur grille complete
            write(nom_fich,999) int(100.*temps),niv_max,s,epsilon
            open(31,file=nom_fich,status='unknown')

            write(31,*) 'TITLE="3D_Dyadic_Tree_OS7DIF2"'
            write(31,*) 'VARIABLES="x","y","z","rho",', &
                '"u_x","u_y","u_z","P","T","w_x","w_y","w_z","helicite"'
            write(31,*) 'ZONE T="Pure_BRICK"', &
                    ', NODES = ', nbre_coins, &
                    ', ELEMENTS = ', compteur_feuille, &
                    ', DATAPACKING=BLOCK, ZONETYPE=FEBRICK', &
                    ', VARLOCATION=([1-3]=NODAL, [4-13]=CELLCENTERED)'
! .....ecriture de la liste des coins des cellules
            write(31, 900)  ((x_coins(l, m), l = 1, nbre_coins), &
                              m = 1, ndim )

! .....ecriture des variables au centre des cellules
            write(31, 900) ((var_print(l, m), l = 1, compteur_feuille), &
                             m = 1, nprint_var)

! .....ecriture de la liste de connectivite
            do l = 1, compteur_feuille
              write(31, 800)  (ind_coins(l, m), m = 1, 2**ndim)
            enddo

            close(31)

          endif

          DEALLOCATE( flag_coins, x_coins, ind_coins, var_print, coords )

        endif

! +++++Algorithme de Strang pour la tensorisation directionnelle

        do n_st = 1, strang_pas

          n_ite = nbre_iter + n_st + (nt - 1) * strang_pas
!          print*, 'n_ite',n_ite,'n_t',nt

! +++++Calcul du pas de temps

! .....initialisation du pas de temps convection et acoustique
          dt = 1.e+10
          dt_ac = 1.e+10

          call pas_de_temps_sp

          print*
          print*, "*****************************************************************"
          write(*,'(a,e13.6,x,a,e13.6,x,a,i3.3)') "dt_ac = ",dt_ac,"dts2 = ",dts2,"m_ac = ", m_ac
		  write(*,'(a,x,e13.6)')"dt = ",dt

! .....accumulation du temps de simulation

          if( (temps+dt) > tps_deb ) dt = abs(tps_deb - temps)

          temps = temps + dt

! .....Mise a jour de la solution (u_ndt) au temps n * Delta t

          do k = 1, nbre_racine(3)
            do j = 1, nbre_racine(2)
              do i = 1, nbre_racine(1)

                root => racine(i, j, k)%ptr
                call undt_update(root)

              enddo
            enddo
          enddo

! +++++integration sur les feuilles uniquement non fictives

          do l = 1, ndim

! .....choix de la direction
            dir = strang(l,n_st)
		  print*, 'dir = ',dir
! +++++Integration Euler + Correction MP

! .....Mise a jour de la solution (u0) pour integration M-C

            do k = 1, nbre_racine(3)
              do j = 1, nbre_racine(2)
                do i = 1, nbre_racine(1)

                  root => racine(i, j, k)%ptr
                  call u0_update(root)

                enddo
              enddo
            enddo

! .....calcul des flux aux interfaces : Euler + Correction OSMP7
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! +++++++++++++++++++++++++++++++++++++++++++++++
! +++++calcul intermédiaire sur dt/2 : convection
! +++++++++++++++++++++++++++++++++++++++++++++++

         call flux_osmp7_conv(dir)

! .....Integration : Euler + Correction

	    call integration_euler(dir,dts2)

! .....Transfert unp1 ==> u

         call mise_a_jour

! .....Codage par valeurs moyennes

            do k = 1, nbre_racine(3)
              do j = 1, nbre_racine(2)
                do i = 1, nbre_racine(1)

                  root => racine(i, j, k)%ptr
                  call encoding_moy(root)

                enddo
              enddo
            enddo

! .....Calcul des valeurs de details sur tout l'arbre
            det_max = 0.

            do k = 1, nbre_racine(3)
              do j = 1, nbre_racine(2)
                do i = 1, nbre_racine(1)

                  root => racine(i, j, k)%ptr
                  call encoding_det(root)

                enddo
              enddo
            enddo

! .....Mise a jour des feuilles fictives

            do k = 1, nbre_racine(3)
              do j = 1, nbre_racine(2)
                do i = 1, nbre_racine(1)

                  root => racine(i, j, k)%ptr
                  call valeur_feuille_fictive(root)

                enddo
              enddo
            enddo

		  do k = 1, nbre_racine(3)
              do j = 1, nbre_racine(2)
                do i = 1, nbre_racine(1)

                  root => racine(i, j, k)%ptr
                  call u0_update(root)

                enddo
              enddo
            enddo

! +++++++++++++++++++++++++++++++++++++++++++++
! +++++calcul intermédiaire sur dt : acoustique
! +++++++++++++++++++++++++++++++++++++++++++++

do i_sp = 1, m_ac

	    call flux_osmp7_acoustic(dir)

! .....Integration : Euler + Correction

	    call integration_euler(dir,dt_ac)

! .....Transfert unp1 ==> u

         call mise_a_jour

! .....Codage par valeurs moyennes

            do k = 1, nbre_racine(3)
              do j = 1, nbre_racine(2)
                do i = 1, nbre_racine(1)

                  root => racine(i, j, k)%ptr
                  call encoding_moy(root)

                enddo
              enddo
            enddo

! .....Calcul des valeurs de details sur tout l'arbre
            det_max = 0.

            do k = 1, nbre_racine(3)
              do j = 1, nbre_racine(2)
                do i = 1, nbre_racine(1)

                  root => racine(i, j, k)%ptr

                  call encoding_det(root)

                enddo
              enddo
            enddo

! .....Mise a jour des feuilles fictives

            do k = 1, nbre_racine(3)
              do j = 1, nbre_racine(2)
                do i = 1, nbre_racine(1)

                  root => racine(i, j, k)%ptr

                  call valeur_feuille_fictive(root)

                enddo
              enddo
            enddo

            do k = 1, nbre_racine(3)
              do j = 1, nbre_racine(2)
                do i = 1, nbre_racine(1)

                  root => racine(i, j, k)%ptr
                  call u0_update(root)

                enddo
              enddo
            enddo

end do

! !
! ! +++++++++++++++++++++++++++++++++++++++++++++++
! ! +++++calcul intermédiaire sur dt/2 : convection
! ! +++++++++++++++++++++++++++++++++++++++++++++++
! !
          call flux_osmp7_conv(dir)

! .....Integration : Euler + Correction

	     call integration_euler(dir,dts2)

! .....Transfert unp1 ==> u

          call mise_a_jour

! .....Codage par valeurs moyennes

            do k = 1, nbre_racine(3)
              do j = 1, nbre_racine(2)
                do i = 1, nbre_racine(1)

                  root => racine(i, j, k)%ptr
                  call encoding_moy(root)

                enddo
              enddo
            enddo

! .....Calcul des valeurs de details sur tout l'arbre
            det_max = 0.

            do k = 1, nbre_racine(3)
              do j = 1, nbre_racine(2)
                do i = 1, nbre_racine(1)

                  root => racine(i, j, k)%ptr

                  call encoding_det(root)

                enddo
              enddo
            enddo

! .....Mise a jour des feuilles fictives

            do k = 1, nbre_racine(3)
              do j = 1, nbre_racine(2)
                do i = 1, nbre_racine(1)

                  root => racine(i, j, k)%ptr

                  call valeur_feuille_fictive(root)

                enddo
              enddo
            enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! +++++ Flux visqueux -  si simulation N-S

            if( Reynolds /= 0. ) then

! +++++Predicteur

! .....Mise a jour de la solution (u0) pour integration Predicteur

              do k = 1, nbre_racine(3)
                do j = 1, nbre_racine(2)
                  do i = 1, nbre_racine(1)

                    root => racine(i, j, k)%ptr
                    call u0_update(root)

                  enddo
                enddo
              enddo

! .....Calcul des flux visqueux aux interfaces

              call flux_visc(dir)

! .....Integration : Predicteur Visqueux

              call integration1_visqueux(dir)

! .....Transfert unp1 ==> u

              call mise_a_jour

! .....Codage par valeurs moyennes

              do k = 1, nbre_racine(3)
                do j = 1, nbre_racine(2)
                  do i = 1, nbre_racine(1)

                    root => racine(i, j, k)%ptr
                    call encoding_moy(root)

                  enddo
                enddo
              enddo

! .....Calcul des valeurs de details sur tout l'arbre
              det_max = 0.

              do k = 1, nbre_racine(3)
                do j = 1, nbre_racine(2)
                  do i = 1, nbre_racine(1)

                    root => racine(i, j, k)%ptr

                    call encoding_det(root)

                  enddo
                enddo
              enddo

! .....Valorisation des feuilles fictives

              do k = 1, nbre_racine(3)
                do j = 1, nbre_racine(2)
                  do i = 1, nbre_racine(1)

                    root => racine(i, j, k)%ptr

                    call valeur_feuille_fictive(root)

                  enddo
                enddo
              enddo

! +++++Correcteur

! .....calcul des flux visqueux aux interfaces

              call flux_visc(dir)

! .....Integration : Correcteur Visqueux

              call integration2_visqueux(dir)

! .....Transfert unp1 ==> u

              call mise_a_jour

! .....Codage par valeurs moyennes

              do k = 1, nbre_racine(3)
                do j = 1, nbre_racine(2)
                  do i = 1, nbre_racine(1)

                    root => racine(i, j, k)%ptr
                    call encoding_moy(root)

                  enddo
                enddo
              enddo

! .....Calcul des valeurs de details sur tout l'arbre
              det_max = 0.

              do k = 1, nbre_racine(3)
                do j = 1, nbre_racine(2)
                  do i = 1, nbre_racine(1)

                    root => racine(i, j, k)%ptr

                    call encoding_det(root)

                  enddo
                enddo
              enddo

! .....Valorisation des feuilles fictives

              do k = 1, nbre_racine(3)
                do j = 1, nbre_racine(2)
                  do i = 1, nbre_racine(1)

                    root => racine(i, j, k)%ptr

                    call valeur_feuille_fictive(root)

                  enddo
                enddo
              enddo

! .....Fin du traitement des flux visqueux
            endif

! .....Fin de boucle sur les directions

! +++++Calcul des residus et valeurs integrales

        Volume = 0.
        errol2(:)=0.
        integral(:)=0.

        ijk_resl2max(:,:) = 0
        err_max(:) = 0.

        call residu

! ! .....ecriture a chaque iteration de direction
! 
!         Print*,' iterations :', n_ite,' ; temps = ',temps
!         write(*,'(a,i6.6,x,a,f8.4,a)') '           Nbre feuilles : ', compteur_feuille, &
!                 ' Compression : ', &
!                 100.*float(compteur_feuille)/float(nbre_maille), ' %'
! 
!         do m = 1, nvar
!           write(*,'(a,i1.1,a,2(x,e15.8),3(x,i4.1))')' Res L2 ( m = ',m,' ) : ', &
!                  sqrt(errol2(m))/float(compteur_feuille), &
!                  err_max(m), ijk_resl2max(:, m)
!         enddo
! 
!         write(*,'(a,5(x,e15.8))')' Int. u dx = ', (integral(m)/Volume, m = 1,nvar)
! 
!         write(15,997) temps, ( integral(m), m = 1, nvar )
! 
! ! +++++ecriture des residus
! 
!         write(10,998) n_ite,temps, &
!              (sqrt(errol2(m))/float(compteur_feuille), m = 1, nvar)


          enddo

! +++++Ecriture et sauvegarde de la solution aux temps intermediaires

          if( (temps >= tps_deb ) .or. &
              (abs(temps-tps_deb) <= zero) ) then

            print*,* 'ecriture dans les fichiers ...'
            Print*,' iterations :', n_ite,' ; temps = ',temps
            write(*,'(a,i6.6,x,a,f8.4,a)') '      Nbre feuilles : ', compteur_feuille, &
                ' Compression : ', &
                 100.*float(compteur_feuille)/float(nbre_maille), ' %'

! .....Allocation des tableaux de travail pour l'ecriture
            ALLOCATE( flag_coins(ijk_max), x_coins(ijk_max,ndim) )

! .....allocation des tableaux coordonnees et champs
            ALLOCATE( ind_coins(compteur_feuille, 2**ndim) )
            ALLOCATE( var_print(compteur_feuille, nprint_var) )

            ALLOCATE( coords(compteur_feuille, ndim) )

! .....Fichier coordonnees des feuilles de l'arbre

            nbre_coins = 0

            flag_coins(:) = -1

            nbre_err = 0
            erreur_linf(:) = 0.
            erreur_l1(:) = 0.
            erreur_l2(:) = 0.

            call print_leaves

            if( ndim == 1 ) then

! .....Fichier des champs sur grille complete
              write(nom_fich,999) int(100.*temps),niv_max,s,epsilon
              open(31,file=nom_fich,status='unknown')

! .....Ecriture au format colonne pour Gnuplot
              write(31,*) '# TITLE="1D_Dyadic_Tree"'
              write(31,*) '# VARIABLES="x","dx","niv","rho","u_x",', &
                          '"P","T"'

              do l = 1, compteur_feuille
                write(31,900) x_coins(l,1), &
                              (var_print(l,m), m = 1, nprint_var)
              enddo

			close(31)

! .....Ecriture au format E-F pour Tecplot ou Paraview
            elseif( ndim == 2 ) then
!
!               write(31,*) 'TITLE="2D_Dyadic_Tree_OS7DIF2"'
!               write(31,*) 'VARIABLES="x","y","rho",', &
!                           '"u_x","u_y","P","T","w"'
!               write(31,*) 'ZONE T="Pure_QUAD"', &
!                     ', NODES = ', nbre_coins, &
!                     ', ELEMENTS = ', compteur_feuille, &
!                     ', DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL', &
!                     ', VARLOCATION=([1-2]=NODAL, [3-8]=CELLCENTERED)'
!
! ! .....ecriture de la liste des coins des cellules
!               write(31, 900)  ((x_coins(l, m), l = 1, nbre_coins), &
!                                 m = 1, ndim )
!
! ! .....ecriture des variables au centre des cellules
!               write(31, 900) ((var_print(l, m), l = 1, compteur_feuille), &
!                                m = 1, nprint_var)
!
! ! .....ecriture de la liste de connectivite
!               do l = 1, compteur_feuille
!                 write(31, 800)  (ind_coins(l, m), m = 1, 2**ndim)
!               enddo

               call save_VTK_2D

! .....Ecriture au format E-F pour Tecplot ou Paraview
            else

! .....Fichier des champs sur grille complete
              write(nom_fich,999) int(100.*temps),niv_max,s,epsilon
              open(31,file=nom_fich,status='unknown')

              write(31,*) 'TITLE="3D_Dyadic_Tree_OS7DIF2"'
              write(31,*) 'VARIABLES="x","y","z","rho",', &
                '"u_x","u_y","u_z","P","T","w_x","w_y","w_z","helicite"'
              write(31,*) 'ZONE T="Pure_BRICK"', &
                    ', NODES = ', nbre_coins, &
                    ', ELEMENTS = ', compteur_feuille, &
                    ', DATAPACKING=BLOCK, ZONETYPE=FEBRICK', &
                    ', VARLOCATION=([1-3]=NODAL, [4-13]=CELLCENTERED)'
! .....ecriture de la liste des coins des cellules
              write(31, 900)  ((x_coins(l, m), l = 1, nbre_coins), &
                                m = 1, ndim )

! .....ecriture des variables au centre des cellules
              write(31, 900) ((var_print(l, m), l = 1, compteur_feuille), &
                               m = 1, nprint_var)

! .....ecriture de la liste de connectivite
              do l = 1, compteur_feuille
                write(31, 800)  (ind_coins(l, m), m = 1, 2**ndim)
              enddo

		    close(31)

            endif

            DEALLOCATE( flag_coins, x_coins, ind_coins, var_print,coords )

            Print*, ' '
            Print*, ' Erreurs :'
            do m = 1, nvar-1
              Print*, ' m = ', m, ' Erreurs Linf = ', erreur_linf(m), &
                      ' Erreurs L1 = ', erreur_l1(m), &
                      ' Erreurs L2 = ', erreur_l2(m)
            enddo

! *****Ecriture de la solution pour sauvegarde

            open(20, file='save_mrchorus.dat', form='unformatted')

            write(20) n_ite, temps, compteur_feuille

            Print*
            Print*, ' Ecriture solution : Iterations = ', n_ite
            Print*, '                          Temps = ', temps
            Print*, '               Nbre de Feuilles = ', compteur_feuille

            call print_reprise

            close(20)

            if( tps_deb >= tps_fin ) cok = .true.

            tps_deb=tps_deb+dtinc

          endif

! +++++Fin de calcul si condition realisee
          if( cok ) goto 69

! .....Fin de boucle sur la tensorisation
        enddo

        print*, "*****************************************************************"
	   print*

! .....Fin de boucle en temps

      Enddo

69    continue


	 call cpu_time(TempsFin)

      print*, 'Temps calcul total : ',TempsFin-TempsDeb

! *****Ecriture de la solution finale (si ce n'est deja fait)

      if( .not.cok ) then

        open(20, file='save_mrchorus.dat', form='unformatted')

        write(20) n_ite, temps, compteur_feuille

        Print*
        Print*, ' Ecriture solution : Iterations = ', n_ite
        Print*, '                          Temps = ', temps
        Print*, '               Nbre de Feuilles = ', compteur_feuille

        call print_reprise

        close(20)

      endif

      close(10)
      close(15)

      open(2,file='N.ech',form='formatted',position='rewind')
	 write(2,'(i3.3)') ech+1
	 close(2)


999   format('Tree_t',i5.5,'_niv',i2.2,'_rhoL1s',i1,'_eps',1pe6.0,'.dat')

998   format(1x,i8,5(1x,1pe15.8))

997   format(5(1x,1pe15.8))

995   format('Int_OS7DIF2_niv',i2.2,'_rhoL1s',i1,'_eps',1pe6.0,'.dat')

900   format(20(1x,1pe22.15))

800   format(8(1x,i6))

      end program MR_CHORUS
