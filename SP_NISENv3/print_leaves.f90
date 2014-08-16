      subroutine print_leaves

      use mod_common
      use mod_structure
      use mod_fonction

      IMPLICIT NONE
   
      INTEGER :: i, j, k, l, leaf, m
      INTEGER :: ijk, Jmj
      INTEGER, DIMENSION(3) :: i_c
      INTEGER, DIMENSION(8,3) :: ordre_coins

      REAL(kind=8) :: rho_ec, PsgM2, TsgM2, helicite
      REAL(kind=8), DIMENSION(3) :: vorticite
      REAL(kind=8), DIMENSION(ndim) :: grad_rho, grad_P, grad_T
      REAL(kind=8), DIMENSION(ndim,ndim) :: grad_u

      REAL(kind=8) :: rayon, ampli, u_theta, t_loc, delta_temp
      REAL(kind=8) :: rho_exact
      REAL(kind=8), DIMENSION(ndim) :: x_0, u_exact

      TYPE(structure_maille), pointer :: maille

! +++++definition de l'ordre de parcours des coins
      ordre_coins(1,1) = 0
      ordre_coins(1,2) = 0
      ordre_coins(1,3) = 0
      ordre_coins(2,1) = 1
      ordre_coins(2,2) = 0
      ordre_coins(2,3) = 0
      ordre_coins(3,1) = 1
      ordre_coins(3,2) = 1
      ordre_coins(3,3) = 0
      ordre_coins(4,1) = 0
      ordre_coins(4,2) = 1
      ordre_coins(4,3) = 0
      ordre_coins(5,1) = 0
      ordre_coins(5,2) = 0
      ordre_coins(5,3) = 1
      ordre_coins(6,1) = 1
      ordre_coins(6,2) = 0
      ordre_coins(6,3) = 1
      ordre_coins(7,1) = 1
      ordre_coins(7,2) = 1
      ordre_coins(7,3) = 1
      ordre_coins(8,1) = 0
      ordre_coins(8,2) = 1
      ordre_coins(8,3) = 1

! +++++Eciture au format colonne pour le 1D

      if( ndim == 1 ) then

! *****On parcourt les feuilles reelles (non-fictives)
        Do leaf = 1, compteur_feuille

          maille => feuille_reelle(leaf)%ptr

          rho_ec = 0.
          do m = 2,nvar-1
            rho_ec = rho_ec + maille%u(m)**2
          enddo
          rho_ec = 0.5 * rho_ec / maille%u(1)

          PsgM2 = (gamma-1.) * ( maille%u(nvar) - rho_ec )

          TsgM2 = PsgM2 / maille%u(1)

! .....Calcul des gradients
          call gradient(maille, grad_rho, grad_u, grad_P, grad_T)

! +++++Ecriture des coordonnees de l'arbre et des variables
          do m = 1, ndim
            x_coins(leaf, m) = maille%x(m)
          enddo

          var_print(leaf, 1) = maille%dx(1)
          var_print(leaf, 2) = maille%maille_niv
          var_print(leaf, 3) = maille%u(1)
          do m = 1, ndim
            var_print(leaf, m+3) = maille%u(m+1)/maille%u(1)
          enddo
          var_print(leaf, ndim+4) = Psgm2
          var_print(leaf, ndim+5) = Tsgm2

        Enddo
    
! +++++Ecriture au format E-F pour le 2-D       

      Elseif( ndim == 2) then

! *****On parcourt les feuilles reelles (non-fictives)

        Do leaf = 1, compteur_feuille

          maille => feuille_reelle(leaf)%ptr
          
          do i = 1, ndim
          coords(leaf,i) = maille%x(i)
          end do

          rho_ec = 0.
          do m = 2,nvar-1
            rho_ec = rho_ec + maille%u(m)**2
          enddo
          rho_ec = 0.5 * rho_ec / maille%u(1)

          PsgM2 = (gamma-1.) * ( maille%u(nvar) - rho_ec )

          TsgM2 = PsgM2 / maille%u(1)

! .....Calcul des gradients
          call gradient(maille, grad_rho, grad_u, grad_P, grad_T)

! .....Calcul de l'erreur : cas particulier du transport du tourbillon

          ampli = 5.

          x_0(1) = temps * v_left(1)
          if( temps * v_left(1) > lim_fin(1) ) &
            x_0(1) = temps * v_left(1) - (lim_fin(1)-lim_deb(1))

          x_0(2) = temps * v_left(2)
          if( temps * v_left(2) > lim_fin(2) ) &
            x_0(2) = temps * v_left(2) - (lim_fin(2)-lim_deb(2))

          rayon = 0.
          do m = 1,ndim
            rayon = rayon + ( maille%x(m)-x_0(m) )**2
          enddo
          rayon = sqrt(rayon)

          delta_temp = -(gamma-1.)*ampli**2*exp(1.-rayon**2) / &
                        (8.*gamma*pi**2) 

          u_theta = 0.5*ampli*exp(.5*(1.-rayon**2)) / pi

          t_loc = 1. + delta_temp * gamma*Mach**2

          rho_exact = rho_left * t_loc**(1./(gamma-1.))

          u_exact(1) = rho_exact * ( v_left(1) - &
                                     u_theta*(maille%x(2)-x_0(2)) )
          u_exact(2) = rho_exact * ( v_left(2) + &
                                     u_theta*(maille%x(1)-x_0(1)) )

          nbre_err = nbre_err + 1

          erreur_linf(1) = max( erreur_linf(1), abs(maille%u(1)-rho_exact) )
          do m = 2, nvar-1
            erreur_linf(m) = max( erreur_linf(m) , &
                           abs(maille%u(m)-rho_exact*u_exact(m-1)) )
          enddo

          erreur_l1(1) = erreur_l1(1) + abs(maille%u(1)-rho_exact)
          do m = 2, nvar-1
            erreur_l1(m) = erreur_l1(m) + &
                           abs(maille%u(m)-rho_exact*u_exact(m-1))
          enddo

          erreur_l2(1) = erreur_l2(1) + (maille%u(1)-rho_exact)**2
          do m = 2, nvar-1
            erreur_l2(m) = erreur_l2(m) + &
                           (maille%u(m)-rho_exact*u_exact(m-1))**2
          enddo

! .....Vorticite
          vorticite(1) = 0.
          vorticite(2) = 0.
          vorticite(3) = grad_u(2,1) - grad_u(1,2)

! *****Ecriture des coordonnees de l'arbre et des variables

          Jmj = niv_max - maille%maille_niv

! *****Pour les 2**ndim coins du volume de controle
          do l = 1, 2**ndim

! .....calcul du mono-indice
            i_c(:) = 0
            do m = 1, ndim
              i_c(m) = 2**Jmj*maille%maille_i(m) + &
                       (ordre_coins(l,m)-1)*2**Jmj
            enddo
            ijk = i_c(1)+1+i_c(2)*(i_max+1)+i_c(3)*(i_max+1)*(j_max+1)

! .....Ecriture des elements
            if( flag_coins(ijk) < 0 ) then
              nbre_coins = nbre_coins + 1
              flag_coins(ijk) = nbre_coins
              do m = 1, ndim
                x_coins(nbre_coins, m) = maille%x(m) + &
                  (real(ordre_coins(l,m))-0.5)*maille%dx(m)
              enddo
            endif
            ind_coins(leaf, l) = flag_coins(ijk)

          enddo

! .....Ecriture des variables
          var_print(leaf, 1) = maille%u(1)
          do m = 1, ndim
            var_print(leaf, m+1) = maille%u(m+1)/maille%u(1)
          enddo
          var_print(leaf, ndim+2) = PsgM2
          var_print(leaf, ndim+3) = TsgM2
          var_print(leaf, ndim+4) = vorticite(3)

        Enddo

! +++++Ecriture au format E-F pour le 3-D       

      Else

! *****On parcourt les feuilles reelles (non-fictives)

        Do leaf = 1, compteur_feuille

          maille => feuille_reelle(leaf)%ptr

          rho_ec = 0.
          do m = 2,nvar-1
            rho_ec = rho_ec + maille%u(m)**2
          enddo
          rho_ec = 0.5 * rho_ec / maille%u(1)

          PsgM2 = (gamma-1.) * ( maille%u(nvar) - rho_ec )

          TsgM2 = PsgM2 / maille%u(1)

! .....Calcul des gradients
          call gradient(maille, grad_rho, grad_u, grad_P, grad_T)

! .....Vorticite
          vorticite(1) = grad_u(3,2) - grad_u(2,3)
          vorticite(2) = grad_u(1,3) - grad_u(3,1)
          vorticite(3) = grad_u(2,1) - grad_u(1,2)

! .....Helicite et amplitude de vorticite (omega)
          helicite = 0.
          do m = 1, ndim
            helicite = helicite + maille%u(m+1)*vorticite(m)/maille%u(1)
          enddo

! *****Ecriture des coordonnees de l'arbre et des variables

          Jmj = niv_max - maille%maille_niv

! *****Pour les 2**ndim coins du volume de controle
          do l = 1, 2**ndim

! .....calcul du mono-indice
            i_c(:) = 0
            do m = 1, ndim
              i_c(m) = 2**Jmj*maille%maille_i(m) + &
                       (ordre_coins(l,m)-1)*2**Jmj
            enddo
            ijk = i_c(1)+1+i_c(2)*(i_max+1)+i_c(3)*(i_max+1)*(j_max+1)

! .....Ecriture des elements
            if( flag_coins(ijk) < 0 ) then
              nbre_coins = nbre_coins + 1
              flag_coins(ijk) = nbre_coins
              do m = 1, ndim
                x_coins(nbre_coins, m) = maille%x(m) + &
                  (real(ordre_coins(l,m))-0.5)*maille%dx(m)
              enddo
            endif
            ind_coins(leaf, l) = flag_coins(ijk)

          enddo

! .....Ecriture des variables
          var_print(leaf, 1) = maille%u(1)
          do m = 1, ndim
            var_print(leaf, m+1) = maille%u(m+1)/maille%u(1)
          enddo
          var_print(leaf, ndim+2) = PsgM2
          var_print(leaf, ndim+3) = TsgM2
          do m = 1, ndim
            var_print(leaf, ndim+3+m) = vorticite(m)
          enddo
          var_print(leaf, 2*ndim+4) = helicite

        Enddo

      Endif

      end subroutine print_leaves
