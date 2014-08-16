	subroutine pas_de_temps_sp

	use mod_common
	use mod_structure
     use mod_fonction

	implicit none

        INTEGER :: l, leaf, m

        REAL(kind=8) :: ec, cson, vp, vp2, mu_0

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

! .....Rayon spectral convection et acoustique
          vp = 0.
          vp2 = 0.
          do m = 2,nvar-1
            vp = max( vp, abs(maille%u(m)/ maille%u(1)))
            vp2 = max (vp2, cson)
          enddo

! .....Vicosite moleculaire et conductivite thermique
          call viscosite(maille, mu_0)

! .....Calcul du pas de temps convectif et acoustique avec critere CFL
          if( Reynolds /= 0. ) then
            do l =1, ndim
              dt = min( dt, maille%dx(l)*CFL/vp, &
                    maille%u(1)*maille%dx(l)**2*Reynolds*Prandtl*.90*CFL/ &
                    (gamma*mu_0) )
              dt_ac = min (dt_ac, maille%dx(l)*CFL/vp2, &
				maille%u(1)*maille%dx(l)**2*Reynolds*Prandtl*.90*CFL/ &
                    (gamma*mu_0) )

            enddo
          else
            do l =1, ndim
              dt = min( dt, maille%dx(l)*CFL/vp )
              dt_ac = min( dt_ac, maille%dx(l)*CFL/vp2 )
            enddo
          endif

	Enddo
#IF ACA
	dts2 = dt
	m_ac = ceiling(0.5*dt/dt_ac)
 	dt_ac = 0.5*dt/m_ac
#ENDIF

#IF CAC
	dts2 = 0.5*dt
	m_ac = ceiling(dt/dt_ac)
 	dt_ac = dt/m_ac
#ENDIF

#if AC .or. CA
	dts2 = dt
	m_ac = ceiling(dt/dt_ac)
 	dt_ac = dt/m_ac
#endif

	end subroutine pas_de_temps_sp
