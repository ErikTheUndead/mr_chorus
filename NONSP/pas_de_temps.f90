	subroutine pas_de_temps

	use mod_common
	use mod_structure
        use mod_fonction

	implicit none

        INTEGER :: l, leaf, m

        REAL(kind=8) :: ec, cson, vp, mu_0

	TYPE(structure_maille), pointer :: maille

! +++++On parcourt les feuilles reelles (non-fictives)

        Do leaf = 1, compteur_feuille

          maille => feuille_reelle(leaf)%ptr

! .....energie cinetique
          ec = 0.
          do m = 2,nvar-1
            ec = ec + 0.5*( maille%u(m)/maille%u(1) )**2
          enddo

! .....celerite du son
          cson = sqrt( gamma * (gamma-1.) * &
                 (maille%u(nvar) / maille%u(1) - ec ))

! .....Rayon spectral
          vp = 0.
          do m = 2,nvar-1
            vp = max( vp, abs(maille%u(m)/ maille%u(1)) + cson )
          enddo

! .....Vicosite moleculaire et conductivite thermique
          call viscosite(maille, mu_0)

! .....Calcule du pas de temps avec critere CFL
          if( Reynolds /= 0. ) then
            do l =1, ndim
              dt = min( dt, maille%dx(l)*CFL/vp, &
                    maille%u(1)*maille%dx(l)**2*Reynolds*Prandtl*.90*CFL/ &
                    (gamma*mu_0) )
            enddo
          else
            do l =1, ndim
              dt = min( dt, maille%dx(l)*CFL/vp )
            enddo
          endif

	Enddo

	end subroutine pas_de_temps
