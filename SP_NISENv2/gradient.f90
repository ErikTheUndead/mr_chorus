      subroutine gradient(maille, grad_rho, grad_u, grad_P, grad_T)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      TYPE(structure_maille), pointer :: maille

      INTEGER :: i_dir, m

      REAL(kind=8) :: ec, T_1, T_0, P_1, P_0

      REAL(kind=8), DIMENSION(ndim) :: grad_rho, grad_P, grad_T
      REAL(kind=8), DIMENSION(ndim,ndim) :: grad_u

      TYPE(structure_maille), pointer :: local

! =================== Calcul des Gradients de rho, U, P, T ================

! *****Dans toutes les directions

      do i_dir = 1, ndim

! +++++Point precedent dans la direction consideree

        if( associated(maille%support(i_dir)%precedent) ) then

          local => maille%support(i_dir)%precedent

          grad_rho(i_dir) = (maille%u(1) - local%u(1)) / maille%dx(i_dir)
          
          do m = 1, ndim
            grad_u(m, i_dir) = ( maille%u(m+1)/maille%u(1) - &
                                local%u(m+1)/local%u(1) ) / maille%dx(i_dir)
          enddo

          ec = 0.
          do m = 1, ndim
            ec = ec + maille%u(m+1)**2 
          enddo
          ec = 0.5 * ec / maille%u(1)**2 

          T_1 = (gamma-1.) * ( maille%u(nvar)/maille%u(1) - ec )

          P_1 = maille%u(1) * T_1

          ec = 0.
          do m = 1, ndim
            ec = ec + local%u(m+1)**2 
          enddo
          ec = 0.5 * ec / local%u(1)**2

          T_0 = (gamma-1.) * ( local%u(nvar)/local%u(1) - ec )

          P_0 = local%u(1) * T_0

          grad_T(i_dir) = (T_1 - T_0) / maille%dx(i_dir)
          grad_P(i_dir) = (P_1 - P_0) / maille%dx(i_dir)

!
! +++++Conditions aux limites (sauf periodiques)

        elseif( maille%maille_i(i_dir) == 1 ) then

          grad_rho(i_dir) = 0.
          
          do m = 1, ndim
            grad_u(m, i_dir) = 0.
          enddo
  
! .....Conditions de symetrie ou de paroi
          if( caltype_deb(i_dir) == -1 ) &
            grad_u(i_dir, i_dir) = 2. * maille%u(i_dir+1)/maille%u(1) &
                                      / maille%dx(i_dir)
          if( caltype_deb(i_dir) == 1 ) then
            do m = 1, ndim
              grad_u(m, i_dir) = 2. * maille%u(m+1)/maille%u(1) &
                                    / maille%dx(i_dir)
            enddo
          endif

          grad_T(i_dir) = 0.
          grad_P(i_dir) = 0.

        else

          Print*, ' Gradient : maille precedente inexistante ; niv = ', &
           maille%maille_niv, ' i = ', maille%maille_i
          Stop

        endif

! +++++Point suivant dans la direction consideree

        if( associated(maille%support(i_dir)%suivant) ) then

          local => maille%support(i_dir)%suivant

          grad_rho(i_dir) = 0.5 * ( grad_rho(i_dir) + &
                (local%u(1) - maille%u(1)) / maille%dx(i_dir) )
          
          do m = 1, ndim
            grad_u(m, i_dir) = 0.5 * ( grad_u(m, i_dir) + &
                     (local%u(m+1)/local%u(1) - &
                      maille%u(m+1)/maille%u(1)) / maille%dx(i_dir) )
          enddo

          ec = 0.
          do m = 1, ndim
            ec = ec + local%u(m+1)**2
          enddo
          ec = 0.5 * ec / local%u(1)**2

          T_1 = (gamma-1.) * ( local%u(nvar)/local%u(1) - ec )

          P_1 = local%u(1) * T_1

          ec = 0.
          do m = 1, ndim
            ec = ec + maille%u(m+1)**2
          enddo
          ec = 0.5 * ec / maille%u(1)**2

          T_0 = (gamma-1.) * ( maille%u(nvar)/maille%u(1) - ec )

          P_0 = maille%u(1) * T_0

          grad_T(i_dir) = 0.5 * ( grad_T(i_dir) + &
                          (T_1 - T_0) / maille%dx(i_dir) )
          grad_P(i_dir) = 0.5 * ( grad_P(i_dir) + &
                          (P_1 - P_0) / maille%dx(i_dir) )

!
! +++++Conditions aux limites (sauf periodiques)

        elseif( maille%maille_i(i_dir) == &
                2**maille%maille_niv * nbre_racine(i_dir) ) then

          grad_rho(i_dir) = 0.5 * grad_rho(i_dir) 
          
          do m = 1, ndim
            grad_u(m, i_dir) = 0.5 * grad_u(m, i_dir)
          enddo
  
! .....Conditions de symetrie ou de paroi
          if( caltype_fin(i_dir) == -1 ) &
            grad_u(i_dir, i_dir) =  grad_u(i_dir, i_dir) + &
                 maille%u(i_dir+1)/maille%u(1) / maille%dx(i_dir)

          if( caltype_fin(i_dir) == 1 ) then
            do m = 1, ndim
              grad_u(m, i_dir) = grad_u(m, i_dir) + &
                  maille%u(m+1)/maille%u(1) / maille%dx(i_dir)
            enddo
          endif

          grad_T(i_dir) = 0.5 * grad_T(i_dir)
          grad_P(i_dir) = 0.5 * grad_P(i_dir)

        else

          Print*, ' Gradient : maille suivante inexistante ; niv = ', &
           maille%maille_niv, ' i = ', maille%maille_i
          Stop

        endif

      enddo

      end subroutine gradient
