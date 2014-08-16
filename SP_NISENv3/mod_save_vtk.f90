MODULE mod_save_vtk

  use mod_common
  use mod_structure

  CONTAINS

  subroutine save_1D

  implicit none

  INTEGER :: m,l

  CHARACTER(len=8) iterationname
  CHARACTER(len=3) echantillonname
  CHARACTER(len=40) filename
  CHARACTER(len=256) :: foldername,full_name,cmd

  write(iterationname,'(I8.8)') int(100.*temps)
  write(echantillonname,'(I3.3)') ech

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Création d'un nouveau dossier pour chaque échantillon
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(filename,*) echantillonname,'_1DSP_NISENv3',iterationname,'.dat'
  write(foldername,*) trim(adjustl(nom_fich)),'.folder'
  full_name = 'WORK_DIR/' // trim(adjustl(foldername)) // '/' // trim(adjustl(filename))
  full_name = trim(full_name)
  foldername = trim(adjustl(foldername))
  cmd = 'mkdir -p WORK_DIR/' // foldername
  cmd = adjustl(trim(cmd))
  call system(cmd)

  open(unit=31, file=full_name, status='unknown')

  ! .....Ecriture au format colonne pour Gnuplot
  write(31,*) '# TITLE="1D_Dyadic_Tree"'
  write(31,*) '# VARIABLES="x","dx","niv","rho","u_x","P","T"'

  do l = 1, compteur_feuille
  write(31,900) x_coins(l,1),(var_print(l,m), m = 1, nprint_var)
  enddo

  close(31)

  900   format(20(1x,f22.15))

! 	999   format(A3,'_1DSP_NISENv1','_niv',i2.2,'_rhoL1s',i1,'_eps',1pe6.0,'/')


  end subroutine save_1D

  subroutine save_VTK_2D

  implicit none

  INTEGER :: m,l

  CHARACTER(len=8) iterationname
  CHARACTER(len=3) echantillonname
  CHARACTER(len=40) filename
  CHARACTER(len=256) :: foldername,full_name,cmd

  write(iterationname,'(I8.8)') int(100.*temps)
  write(echantillonname,'(I3.3)') ech

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Création d'un nouveau dossier pour chaque échantillon
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   call system('DIR_WORK=~/duhem/stage_2014_Limsi/WORK_DIR')
  write(filename,*) echantillonname,'_2DSP_NISENv3_',iterationname,'.vtk'
!   write(foldername,*) echantillonname,'_2DSP_NISENv1.folder'
  write(foldername,*) trim(adjustl(nom_fich)),'.folder'
  full_name = 'WORK_DIR/' // trim(adjustl(foldername)) // '/' // trim(adjustl(filename))
  full_name = trim(full_name)
  foldername = trim(adjustl(foldername))
  cmd = 'mkdir -p WORK_DIR/' // foldername
  cmd = adjustl(trim(cmd))
  call system(cmd)
!   call system('DIR_PROG=/people/duhem/stage_2014_Limsi/SP_MR_CHORUS_ISENv1')
!   call system('mkdir -p WORK_DIR')
!   call system('DIR_WORK=WORK_DIR/')
!   write(cmd,*) 'mkdir -p $DIR_WORK/' // trim(adjustl(foldername))
!   cmd = trim(cmd)
!   print*, cmd
!   call system(cmd)
!   call system('mkdir -p $DIR_WORK/' // trim(adjustl(foldername)))

  open(unit=20, file=full_name, status='unknown')

  write(20,'(a)')'# vtk DataFile Version 3.1'
  write(20,'(a)')'2D MR_CHORUS_SPLITTING'
  write(20,'(a)')'ASCII'
  write(20,'(a)')'DATASET UNSTRUCTURED_GRID'
  write(20,*)
! .....ecriture des coordonnées des points
  write(20,'((A7),(I7),(A6))') 'POINTS ',nbre_coins," float"
  do l = 1, nbre_coins
   write(20, 802) (x_coins(l, m), m = 1, ndim), ' 0.0'
  end do
  write(20,*)

  ! .....ecriture des variables
  write(20,*) 'CELLS', compteur_feuille, 5*compteur_feuille
  do l = 1, compteur_feuille
   write(20, 902) '4', (ind_coins(l, m)-1, m = 1, 2**ndim)
  enddo
  write(20,*)

  write(20,*) 'CELL_TYPES', compteur_feuille
  do l = 1, compteur_feuille
  write(20,*) '9'
  end do
  write(20,*)

  write(20,*) 'CELL_DATA', compteur_feuille
  write(20,*) 'SCALARS rho float'
  write(20,*) 'LOOKUP_TABLE default'
  do l = 1, compteur_feuille
  write(20,'(f22.15)') var_print(l,1)
  end do
  write(20,*)

  write(20,*) 'SCALARS pressure float'
  write(20,*) 'LOOKUP_TABLE default'
  do l = 1, compteur_feuille
  write(20,'(f22.15)') var_print(l,4)
  end do
  write(20,*)

  write(20,*) 'SCALARS temperature float'
  write(20,*) 'LOOKUP_TABLE default'
  do l = 1, compteur_feuille
  write(20,'(f22.15)') var_print(l,5)
  end do
  write(20,*)

  write(20,*) 'SCALARS vorticity float'
  write(20,*) 'LOOKUP_TABLE default'
  do l = 1, compteur_feuille
  write(20,'(f22.15)') var_print(l,6)
  end do
  write(20,*)


  write(20,*) 'VECTORS velocity float'
  do l = 1, compteur_feuille
  write(20,'(3(f22.15),X)') var_print(l,2), var_print(l,3), 0.
  end do
  write(20,*)

  close(20)

!   cmd = 'mv ' // trim(adjustl(filename)) // ' $DIR_WORK/' // trim(adjustl(foldername))
!   cmd = trim(cmd)
!   print*, cmd
!   call system(cmd)
! !   call system('mkdir -p WORK_DIR/toto')
! !   mv Echantillon2DSP_ISEN_num_030.folder/ ~/stage_2014_Limsi/WORK_DIR


  902   format(a,4(i6))
  802   format(2(f15.8,1x),a)

  end subroutine save_VTK_2D

end module mod_save_vtk