      subroutine calfin_flux_visc(maille, dir, fvisc_droit)

      use mod_common
      use mod_structure
      use mod_fonction

      IMPLICIT NONE

      TYPE(structure_maille), pointer :: maille

      INTEGER :: dir

      REAL(kind=8), DIMENSION(nvar) :: fvisc_droit

      INTEGER :: l, m

      REAL(kind=8), DIMENSION(nvar, 3) :: u_n

      REAL(kind=8) :: rho_ec
      REAL(kind=8), DIMENSION(3) :: mu, tsgM2
      INTEGER, DIMENSION(4,3) :: ijk_dt

      REAL(kind=8) :: mu_wall
      REAL(kind=8) :: divV, T_xyz
      REAL(kind=8), DIMENSION(ndim, ndim) :: u_xyz

      REAL(kind=8), DIMENSION(3) :: normale

      TYPE(structure_maille), pointer :: courant

! .....definition des indices des cellules a droite :
!      ijk_dt( 2**(ndim-1) , dir )
  
      ijk_dt(1,1) = 2
      ijk_dt(2,1) = 4
      ijk_dt(3,1) = 6
      ijk_dt(4,1) = 8
  
      ijk_dt(1,2) = 3
      ijk_dt(2,2) = 4
      ijk_dt(3,2) = 7
      ijk_dt(4,2) = 8
  
      ijk_dt(1,3) = 5
      ijk_dt(2,3) = 6
      ijk_dt(3,3) = 7
      ijk_dt(4,3) = 8

! ======== Flux Visqueux : Condition a la limite amont du domaine ========

! +++++condition de frontiere libre

      if( caltype_fin(dir) == 0 ) then

        do m = 1, nvar
          fvisc_droit(m) = maille%fvisc_gch(m, dir)
        enddo
        
! +++++condition de paroi

      elseif( caltype_fin(dir) == 1 ) then

! *****Extrapolation de la temperature et viscosite (ordre 2)

! .....1er point au dessus de la paroi
        courant => maille

        u_n(1,1) = courant%u(1)
        rho_ec = 0.
        do m = 2,nvar-1
          u_n(m,1) = courant%u(m)
          rho_ec = rho_ec + courant%u(m)**2
        enddo
        rho_ec = 0.5 * rho_ec / courant%u(1)
        u_n(nvar,1) = courant%u(nvar)

        tsgM2(1) = (gamma-1.)*(courant%u(nvar)-rho_ec) / courant%u(1)

        call viscosite(courant, mu(1))
        
! .....2eme point au dessus de la paroi

        courant => courant%support(dir)%precedent

        u_n(1,2) = courant%u(1)
        rho_ec = 0.
        do m = 2,nvar-1
          u_n(m,2) = courant%u(m)
          rho_ec = rho_ec + courant%u(m)**2
        enddo
        rho_ec = 0.5 * rho_ec / courant%u(1)
        u_n(nvar,2) = courant%u(nvar)

        tsgM2(2) = (gamma-1.)*(courant%u(nvar)-rho_ec) / courant%u(1)

        call viscosite(courant, mu(2))

! .....3eme point au dessus de la paroi
        courant => courant%support(dir)%precedent

        u_n(1,3) = courant%u(1)
        rho_ec = 0.
        do m = 2,nvar-1
          u_n(m,3) = courant%u(m)
          rho_ec = rho_ec + courant%u(m)**2
        enddo
        rho_ec = 0.5 * rho_ec / courant%u(1)
        u_n(nvar,3) = courant%u(nvar)

        tsgM2(3) = (gamma-1.)*(courant%u(nvar)-rho_ec) / courant%u(1)

        call viscosite(courant, mu(3))

! *****Gradient de vitesses
        do l = 1, ndim
          do m = 2, nvar-1
            u_xyz(m-1, l) = 0.
          enddo            
        enddo

! .....On suppose nulle la vitesse a la paroi
        do m = 2, nvar-1
! .....Ordre 2
!!          u_xyz(m-1, dir) = (-9.*u_n(m, 1) / u_n(1, 1) + &
!!                                 u_n(m, 2) / u_n(1, 2))/(3.*maille%dx(dir))
! .....Ordre 3
          u_xyz(m-1, dir) = (-225.*u_n(m, 1) / u_n(1, 1) + &
                               50.*u_n(m, 2) / u_n(1, 2) - &
                                9.*u_n(m, 3) / u_n(1, 3) )/(60.*maille%dx(dir))
        enddo            

! *****Divergence de vitesse
        divV = 0.
        do m = 2, nvar-1
          divV = divV + u_xyz(m-1,m-1) 
        enddo

! ======== Flux Visqueux : Condition a la limite amont du domaine =====

! *****Maille courante

        courant => maille

! .....vicosite dynamique a la paroi

        mu_wall = ( 15.*mu(1) - 10.*mu(2) + 3.*mu(3) )/8.      
          
! .....Composantes du Flux visqueux
        fvisc_droit(1) = 0.
        do m = 2, nvar-1
          fvisc_droit(m) = - mu_wall*( u_xyz(dir,m-1) + u_xyz(m-1,dir) - &
                             2.*uns3*divV*delta_ij(dir,m-1) ) / Reynolds
        enddo 

! .....Paroi adiabatique
        fvisc_droit(nvar) = 0.

! .....Sinon temperature de paroi imposee
        if( .not.cle_twall ) then
 
! .....gradient de Temperature suivant "dir"
        
          T_xyz = ( 8.*t_wall - 9.*tsgM2(1)*gamma*Mach**2 + &
                       tsgM2(2)*gamma*Mach**2 )/(3.*maille%dx(dir))

          fvisc_droit(nvar) = - mu_wall*T_xyz / &
                                ((gamma-1.)*Prandtl*Reynolds*Mach**2)

        endif

      endif

      end subroutine calfin_flux_visc
