      recursive subroutine init_sol(maille)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      INTEGER :: l, m
      REAL(kind=8) :: rho_ec, gh_inf
      REAL(kind=8) :: rhot_inf, pt_inf, p_loc, s_choc
      REAL(kind=8) :: R_0, rayon, ampli, c2, u_theta, t_loc, delta_temp
      REAL(kind=8) :: rho_aval, p_aval
      REAL(kind=8), DIMENSION(ndim) :: x_0, x_low, x_high

      TYPE(structure_maille), pointer:: maille

! *****traitement des feuilles

      if(.not.associated(maille%fils(1)%ptr)) then

! +++++Initialisation globale

        maille%u(1) = rho_left

        rho_ec = 0.
        do m = 2,nvar-1
          maille%u(m) = rho_left*v_left(m-1)
          rho_ec = rho_ec + 0.5*rho_left*v_left(m-1)**2
        enddo

        maille%u(nvar) = p_left/(gamma*(gamma-1.)*Mach**2) + rho_ec

        gh_inf = p_left/((gamma-1.)*rho_left*Mach**2) + 0.5

        rhot_inf = rho_left * &
                   (1.+ 0.5*(gamma-1.)*Mach**2 ) ** (1./(gamma-1.))

        pt_inf = p_left * &
                 (1.+ 0.5*(gamma-1.)*Mach**2 ) ** (gamma/(gamma-1.))

#if Cav1
         if( maille%maille_i(2) == 2**maille%maille_niv*nbre_racine(2) ) then

          u_cav(maille%maille_i(1),1 ) =  4. * (maille%x(1))**2 * (lim_fin(1) - (maille%x(1))**2)
          u_cav(maille%maille_i(1),2 ) = v_right(2)

          maille%u(2) = rho_right*u_cav( maille%maille_i(1), 1 )
          maille%u(3) = rho_right*u_cav( maille%maille_i(1), 2 )
          rho_ec = 0.5*(maille%u(2)**2/rho_right + maille%u(3)**2/rho_right)
          maille%u(nvar) = p_left/(gamma*(gamma-1.)*Mach**2) + rho_ec

        end if
#endif

#if Cav2
         if( maille%maille_i(2) == 2**maille%maille_niv*nbre_racine(2) ) then

          u_cav(maille%maille_i(1),1 ) = v_right(1)
          u_cav(maille%maille_i(1),2 ) = v_right(2)

          maille%u(2) = rho_right*u_cav( maille%maille_i(1), 1 )
          maille%u(3) = rho_right*u_cav( maille%maille_i(1), 2 )
          rho_ec = 0.5*(maille%u(2)**2/rho_right + maille%u(3)**2/rho_right)
          maille%u(nvar) = p_left/(gamma*(gamma-1.)*Mach**2) + rho_ec

        end if
#endif

! +++++initialisation convection tourbillon

#if CT

       ampli = 5.
       x_0(1) = 0.
       x_0(2) = 0.

       rayon = 0.
       do m = 1,ndim
         rayon = rayon + ( maille%x(m)-x_0(m) )**2
       enddo
       rayon = sqrt(rayon)

       delta_temp = -(gamma-1.)*ampli**2*exp(1.-rayon**2) / &
                     (8.*pi**2)

       u_theta = 0.5*ampli*exp(.5*(1.-rayon**2)) / pi
!      t_loc = 1. + delta_temp
       t_loc = 1. + delta_temp * gamma*Mach**2

!      maille%u(1) = p_left/t_loc

       maille%u(1) = rho_left

       maille%u(2) = maille%u(1) * ( v_left(1) - &
                                     u_theta*(maille%x(2)-x_0(2)) )
       maille%u(3) = maille%u(1) * ( v_left(2) + &
                                     u_theta*(maille%x(1)-x_0(1)) )

       rho_ec = 0.
       do m = 2,nvar-1
         rho_ec = rho_ec + maille%u(m)**2
       enddo
       rho_ec = 0.5 * rho_ec / maille%u(1)

       maille%u(nvar) = maille%u(1)*t_loc/(gamma*(gamma-1.)*Mach**2) + &
                        rho_ec
!     maille%u(nvar) = p_left/(gamma*(gamma-1.)*Mach**2) + rho_ec

#endif

#if SDT

! +++++initialisation d'un spot de temperature

       x_0(1) = 0.5
       x_0(2) = 0.5

       rayon = 0.
       do m = 1,ndim
         rayon = rayon + ( maille%x(m)-x_0(m) )**2
       enddo
       rayon = sqrt(rayon)

       c2 = 7.
!        ampli = (rayon**2 - c2**2)**2 / c2**4
       ampli = (rayon**2 - c2**2)**2 / c2**4

       c2 = 1./(0.07**2)

       delta_temp = (gamma -1.)*ampli*exp(1-c2*rayon**2)

!        t_loc = 1. + delta_temp
       t_loc = 1. + delta_temp * gamma*Mach**2

!        if(maille%x(1) < (lim_fin(1)+lim_deb(1))*.5) then

!          maille%u(1) = rho_left/t_loc
         maille%u(1) = rho_left

         do m = 2,nvar-1
           maille%u(m) = maille%u(1) * v_left(m-1)
         enddo

         rho_ec = 0.
         do m = 2,nvar-1
           rho_ec = rho_ec + maille%u(m)**2
         enddo
         rho_ec = 0.5 * rho_ec / maille%u(1)

         maille%u(nvar) = maille%u(1) * t_loc/(gamma*(gamma-1.)*Mach**2) + &
                          rho_ec

!        endif

#endif

! +++++initialisation convection tourbillon isentropique

#if CTI

       ampli = 5.
       x_0(1) = 0.
       x_0(2) = 0.

       rayon = 0.
       do m = 1,ndim
         rayon = rayon + ( maille%x(m)-x_0(m) )**2
       enddo
       rayon = sqrt(rayon)

       delta_temp = -(gamma-1.)*ampli**2*exp(1.-rayon**2) / &
                     (8.*pi**2)

       u_theta = 0.5*ampli*exp(.5*(1.-rayon**2)) / pi
       t_loc = 1. + delta_temp * gamma*Mach**2

       maille%u(1) = rho_left * t_loc**(1./(gamma-1.))

       maille%u(2) = maille%u(1) * ( v_left(1) - &
                                     u_theta*(maille%x(2)-x_0(2)) )
       maille%u(3) = maille%u(1) * ( v_left(2) + &
                                     u_theta*(maille%x(1)-x_0(1)) )

       rho_ec = 0.
       do m = 2,nvar-1
         rho_ec = rho_ec + maille%u(m)**2
       enddo
       rho_ec = 0.5 * rho_ec / maille%u(1)

       maille%u(nvar) = maille%u(1)*t_loc/(gamma*(gamma-1.)*Mach**2) + &
                        rho_ec

#endif

! +++++initialisation tube a choc

#if TAC

        if(maille%x(1) <= (lim_fin(1)+lim_deb(1))*.5) then
          maille%u(1) = rho_left
          rho_ec = 0.
          do m = 2,nvar-1
            maille%u(m) = rho_left*v_left(m-1)
            rho_ec = rho_ec + maille%u(m)**2
          enddo
           rho_ec = 0.5 * rho_ec / maille%u(1)
!		  rho_ec = 0.5 * rho_ec * maille%u(1)
          maille%u(nvar) = p_left/(gamma*(gamma-1.)*Mach**2) + rho_ec
        else
          maille%u(1) = rho_right
          rho_ec = 0.
          do m = 2,nvar-1
            maille%u(m) = rho_right*v_right(m-1)
            rho_ec = rho_ec + maille%u(m)**2
          enddo
           rho_ec = 0.5 * rho_ec / maille%u(1)
!		  rho_ec = 0.5 * rho_ec * maille%u(1)
          maille%u(nvar) = p_right/(gamma*(gamma-1.)*Mach**2) + rho_ec
        endif

#endif

#if ODC

! +++++initialisation onde de choc

       rho_left = rho_right * (gamma+1.)*Mach**2 / &
                              ( (gamma-1.)*Mach**2 + 2. )

       p_left = p_right * (2.*gamma*Mach**2-(gamma-1.)) / (gamma+1.)

       gh_inf = p_right/((gamma-1.)*rho_right*Mach**2)
       do m = 2, nvar-1
         gh_inf = gh_inf + 0.5*v_right(m-1)**2
       enddo

       v_left(1) = rho_right*( v_right(1)- 1. )/rho_left + 1.
       v_left(2) = v_right(2)

       if(maille%x(1) <= 0.1 ) then

         maille%u(1) = rho_left

         rho_ec = 0.
         do m = 2,nvar-1
           maille%u(m) = rho_left*v_left(m-1)
           rho_ec = rho_ec + maille%u(m)**2
         enddo
         rho_ec = 0.5 * rho_ec / maille%u(1)

         maille%u(nvar) = p_left/(gamma*(gamma-1.)*Mach**2) + rho_ec

       else

         maille%u(1) = rho_right

         rho_ec = 0.
         do m = 2,nvar-1
           maille%u(m) = rho_right*v_right(m-1)
           rho_ec = rho_ec + maille%u(m)**2
         enddo
         rho_ec = 0.5 * rho_ec / maille%u(1)

         maille%u(nvar) = p_right/(gamma*(gamma-1.)*Mach**2) + rho_ec

       endif

#endif



#if ODCS

! +++++initialisation onde de choc stationaire

! .....valeur de saut Rankine-Hugoniot

       rho_aval = rho_left * (gamma+1.)*Mach**2 &
                           / ( (gamma-1.)*Mach**2 + 2. )

       p_aval = p_left * ( 2.*gamma*Mach**2 - (gamma-1.) ) / (gamma+1.)

       if(maille%x(1) <= (lim_fin(1)+lim_deb(1))*.5) then

         maille%u(1) = rho_left

         rho_ec = 0.
         do m = 2,nvar-1
           maille%u(m) = rho_left*v_left(m-1)
           rho_ec = rho_ec + maille%u(m)**2
         enddo
         rho_ec = 0.5 * rho_ec / maille%u(1)
         maille%u(nvar) = p_left/(gamma*(gamma-1.)*Mach**2) + rho_ec

       else

         maille%u(1) = rho_aval
         do m = 2,nvar-1
           maille%u(m) = rho_left*v_left(m-1)
         enddo
         maille%u(nvar) = rho_aval * gh_inf - p_aval/(gamma*Mach**2)

       endif

#endif



#if TDT

! +++++tourbillon de Taylor

       x_0(1) = 0.
       x_0(2) = 0.
!!        x_0(1) = 0.5
!!        x_0(2) = 0.5
!!!!        x_0(2) = 1.
!!
       rayon = 0.075
       ampli = .3 * exp(.5) / rayon
!!        ampli = .25 * exp(.5) / rayon
       c2 = 1./(2.*rayon*rayon)

       rayon = 0.
       do m = 1, ndim
         rayon = rayon + ( maille%x(m)-x_0(m) )**2
       enddo
       rayon = sqrt(rayon)

       u_theta = ampli*rayon*exp(-c2*rayon**2)

!!        if(maille%x(1) < (lim_fin(1)+lim_deb(1))*.5) then

         maille%u(1) = rho_left

         maille%u(2) = maille%u(2) - &
                       maille%u(1)*u_theta*(maille%x(2)-x_0(2))/rayon
         maille%u(3) = maille%u(3) + &
                       maille%u(1)*u_theta*(maille%x(1)-x_0(1))/rayon

         rho_ec = 0.
         do m = 2,nvar-1
           rho_ec = rho_ec + maille%u(m)**2
         enddo
         rho_ec = 0.5 * rho_ec / maille%u(1)

         maille%u(nvar) = p_left/(gamma*(gamma-1.)*Mach**2) + rho_ec

!!        endif

#endif

#if PTV

! +++++Plane Taylor Vortex
			maille%u(1) = rho_left

			maille%u(2) = sin(pi*maille%x(1))*cos(pi*maille%x(2))

			maille%u(3) = - cos(pi*maille%x(1))*sin(pi*maille%x(2))

			rho_ec = 0.
         do m = 2,nvar-1
           rho_ec = rho_ec + maille%u(m)**2
         enddo
         rho_ec = 0.5 * rho_ec / maille%u(1)

         p_loc = p_left + 0.25*(cos(2*pi*maille%x(1)) + cos(2*pi*maille%x(2)))*gamma*Mach**2

			maille%u(nvar) = p_loc/(gamma*(gamma-1.)*Mach**2) + rho_ec

#endif

#if ICT1D

! +++++initialisation de l'interaction choc-turbulence 1D

       if(maille%x(1) >= 1.) then
!!        if(maille%x(1) >= -0.8) then
         maille%u(1) = rho_right + 0.2 * sin(5.*maille%x(1))
!!          maille%u(1) = rho_right + 0.2 * sin(5.*pi*maille%x(1))
         rho_ec = 0.
         do m = 2,nvar-1
           maille%u(m) = rho_right*v_right(m-1)
           rho_ec = rho_ec + maille%u(m)**2
         enddo
         rho_ec = 0.5 * rho_ec / maille%u(1)
         maille%u(nvar) = p_right/(gamma*(gamma-1.)*Mach**2) + rho_ec
       else
         maille%u(1) = rho_left
         rho_ec = 0.
         do m = 2,nvar-1
           maille%u(m) = rho_left*v_left(m-1)
           rho_ec = rho_ec + maille%u(m)**2
         enddo
         rho_ec = 0.5 * rho_ec / maille%u(1)
         maille%u(nvar) = p_left/(gamma*(gamma-1.)*Mach**2) + rho_ec
       endif

#endif

      else

! +++++traitement des fils

        if( associated(maille%fils(1)%ptr) ) then
          do l = 1, 2**ndim
            call init_sol(maille%fils(l)%ptr)
          enddo
        endif

      endif

      end subroutine init_sol

