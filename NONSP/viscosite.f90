      subroutine viscosite(maille, mu_loc)

      use mod_common
      use mod_structure

      IMPLICIT NONE

      TYPE(structure_maille), pointer :: maille

      REAL(kind=8), INTENT(OUT) :: mu_loc

      INTEGER :: l, m

      REAL(kind=8) :: ec, T_loc

! +++++Vicosite dynamique

! .....Loi de Sutherland
!!#ifdef VISC_SUTHER

!!      ec = 0.
!!      do m = 2,nvar-1
!!        ec = ec + 0.5 * ( maille%u(m)/maille%u(1) )**2
!!      enddo
!!
!!      T_loc = gamma*(gamma-1.)*Mach**2 * &
!!              ( maille%u(nvar) / maille%u(1) - ec )
!!
!!      mu_loc =  T_loc * sqrt(T_loc) * ( 1. + C_suth ) / &
!!                                      ( T_loc + C_suth )

! .....Loi de viscosite constante
!!#else
!!#ifdef VISC_CST

      mu_loc = 1.

!!#endif
!!#endif

      end subroutine viscosite
