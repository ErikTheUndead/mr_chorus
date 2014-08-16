      subroutine lect_don

      use mod_common

      IMPLICIT NONE

      INTEGER :: l, lect

! +++++lecture des donnees dans le fichier don_euler.dat, unite=lect

      lect=1

#if Nbdim .EQ. 1 .and. SDT
		open(lect,file='don_mrchorus1DSDT.dat',form='formatted')
		write(Cas,'(a3)') 'SDT'
#elif Nbdim .EQ. 1
		open(lect,file='don_mrchorus1D.dat',form='formatted')
      write(Cas,'(a7)') 'default'
#endif

#if Nbdim .EQ. 2 .and. CTI
		open(lect,file='don_mrchorus2DCTI.dat',form='formatted')
		write(Cas,'(a3)') 'CTI'
#elif Nbdim .EQ. 2 .and. CT
		open(lect,file='don_mrchorus2DCT.dat',form='formatted')
		write(Cas,'(a2)') 'CT'
#elif Nbdim .EQ. 2 .and. (Cav1 .or. Cav2)
      open(lect,file='don_mrchorus2DCav.dat',form='formatted')
      write(Cas,'(a3)') 'Cav'
#elif Nbdim .EQ. 2 .and. SDT
		open(lect,file='don_mrchorus2DSDT.dat',form='formatted')
		write(Cas,'(a3)') 'SDT'
#elif Nbdim .EQ. 2 .and. TDT
		open(lect,file='don_mrchorus2DTDT.dat',form='formatted')
		write(Cas,'(a3)') 'TDT'
#elif Nbdim .EQ. 2 .and. PTV
		open(lect,file='don_mrchorus2DPTV.dat',form='formatted')
		write(Cas,'(a3)') 'PTV'
#elif Nbdim .EQ. 2
		open(lect,file='don_mrchorus2D.dat',form='formatted')
		write(Cas,'(a7)') 'default'
#endif

#if Nbdim .EQ. 3
      open(lect,file='don_mrchorus3D.dat',form='formatted')
		write(Cas,'(a7)') 'default'
#endif
!
! .....limites du domaine, nombre de racines
!      et flags des cond. aux limites

      nbre_racine(:) = 1

      do l = 1, ndim
        read(lect,*) lim_deb(l), lim_fin(l), nbre_racine(l)
        read(lect,*) caltype_deb(l), caltype_fin(l)
      enddo

! .....initialisation de la solution a gauche
      read(lect,*) rho_left

      do l = 2, nvar-1
        read(lect,*) v_left(l-1)
      enddo

      read(lect,*) p_left

! .....initialisation de la solution a droite
      read(lect,*) rho_right

      do l = 2, nvar-1
        read(lect,*) v_right(l-1)
      enddo

      read(lect,*) p_right

!!
!       p_left = gamma*p_left
!       p_right = gamma*p_right
!!

! .....Longueur de reference
      read(lect,*) L_ref

      if( L_ref == 0. ) L_ref = 1.

! .....Cle condition sur temperature de paroi
      read(lect,*) cle_twall

! .....Temperature de paroi
      read(lect,*) t_wall

! .....Nombre de Mach
      read(lect,*) Mach

! .....Nombre de Reynolds
      read(lect,*) Reynolds

! .....Nombre de Prandtl
      read(lect,*) Prandtl

! .....nombre CFL
      read(lect,*) CFL

! .....nombre d'iteration en temps
      read(lect,*) npdt

! .....choix du nombre de niveaux de grilles
      read(lect,*) niv_max

! .....choix du niveaux minimum de grilles
      read(lect,*) niv_min

! .....choix de l'ordre de reconstruction s
      read(lect,*) s

! .....choix de la valeur d'epsilon
      read(lect,*) epsilon

! .....Cle de reprise ==> reprise = True
      read(lect,*) cle_reprise

! .....premier temps de sauvegarde
      read(lect,*) tps_deb

! .....increment en temps pour la sauvegarde
      read(lect,*) dtinc

! .....temps final de simulation
      read(lect,*) tps_fin

      close(lect)

! +++++Valeurs de reference
		write(*,'(a,i3.3)') 'echantillon',ech

		write(nom_fich,1000) ech,Ndim,trim(Cas),Mach,epsilon,CFL,niv_max,tps_fin
		print*, adjustl(trim(nom_fich))

		open(15, file=nom_fich,form='formatted')!

		write(15,*)
      write(15,*) ' nombre de racines = ', nbre_racine


      if( Reynolds == 0. ) then
        write(15,*) ' Calcul Euler : ', ndim, ' D'
        write(15,*) ' ------------ '
        write(15,*)
        write(15,*) ' - Nombre de Mach      = ', Mach
      else
        write(15,*)
        write(15,*)
        write(15,*) ' Calcul Navier-Stokes : ', ndim, ' D'
        write(15,*), ' -------------------- '
        write(15,*)
        write(15,*) ' - Nombre de Mach      = ', Mach
        write(15,*) ' - Nombre de Reynolds  = ', Reynolds
        write(15,*) ' - Nombre de Prandtl   = ', Prandtl
      endif

! .....nombre total de mailles

      nbre_maille = 2**(ndim * niv_max)
      Nbre_pt_arbre = ( 2**((niv_max+1)*ndim) - 1 ) / ( 2**ndim - 1 )

      do l = 1, ndim
        nbre_maille = nbre_maille * nbre_racine(l)
        Nbre_pt_arbre = Nbre_pt_arbre * nbre_racine(l)
      enddo

      write(15,*)
      write(15,*) ' Nombre total de points dans l''arbre = ', Nbre_pt_arbre
      write(15,*) ' Nombre total de feuilles possibles     = ', nbre_maille

		close(15)

1000	format(i3.3,'_',i1.1,'DNISENv2-',a,'_Mach=',1pe6.0,'_eps=',1pe6.0,'_CFL=',1pe6.0,'_nivmax=',i2.2,'_tmax=',1pe6.0,'.out')

      end subroutine lect_don