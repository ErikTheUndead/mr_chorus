!	module de declaration des variables communes

      module mod_common

      implicit none

      INTEGER, PARAMETER :: ndim = Nbdim
		INTEGER, PARAMETER :: nvar = Nbdim+2
		CHARACTER(len=10) :: Cas

      INTEGER, PARAMETER :: nordre = 4
      INTEGER, PARAMETER :: niveau_maximum = 12

      REAL(kind=8), PARAMETER :: zero = 1.e-14
      REAL(kind=8), PARAMETER :: gamma = 1.4
      REAL(kind=8), PARAMETER :: Rgaz = 287.
      REAL(kind=8), PARAMETER :: cv = Rgaz / (gamma-1.)
      REAL(kind=8), PARAMETER :: cp = gamma * cv
      REAL(kind=8), PARAMETER :: unstref = 1. / 273.
      REAL(kind=8), PARAMETER :: C_suth = 110.4 * unstref
      REAL(kind=8), PARAMETER :: uns3 = 1./3.

      REAL(kind=8), dimension(ndim) :: centre_0
      REAL(kind=8) :: rayon_0

      INTEGER :: niv_max, niv_min, npdt, nbre_iter
      INTEGER :: nt, s, sgrad
      INTEGER :: nbre_maille, Nbre_pt_arbre
      INTEGER :: compteur_create, compteur_feuille, compteur_fictive
      INTEGER, DIMENSION(ndim) :: caltype_deb, caltype_fin
      INTEGER, DIMENSION(3) :: nbre_racine

      INTEGER, DIMENSION(ndim, nvar) :: ijk_resl2max

      REAL(kind=8), DIMENSION(ndim, ndim) :: delta_ij

      REAL(kind=8) :: rho_left, rho_right, p_left, p_right, t_wall
      REAL(kind=8), DIMENSION(ndim) :: v_left, v_right
      REAL(kind=8) :: L_ref, Mach, Reynolds, Prandtl
      REAL(kind=8), DIMENSION(ndim) :: lim_deb, lim_fin, dxmin
      REAL(kind=8) :: temps, dt, tps_deb, dtinc, tps_fin
      REAL(kind=8) :: Volume, det_max, epsilon, pi, CFL
      REAL(kind=8), DIMENSION(nvar) :: errol2, integral

      REAL(kind=8), DIMENSION(nvar) :: err_max

      REAL(kind=8), DIMENSION(:), POINTER :: coef

      REAL(kind=8), DIMENSION(nvar) :: erreur_linf, erreur_l1, erreur_l2
      INTEGER :: nbre_err

! .....definition des variables pour ecriture
      INTEGER :: nbre_coins, i_max, j_max, k_max, ijk_max, nprint_var
      INTEGER, allocatable, dimension(:) :: flag_coins
      INTEGER, allocatable, dimension(:,:) :: ind_coins

      REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: x_coins
      REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: var_print
		CHARACTER(len=80) :: nom_fich

      LOGICAL :: cle_reprise
      LOGICAL :: cle_twall

      Integer :: ech

      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: u_cav

      end module mod_common
