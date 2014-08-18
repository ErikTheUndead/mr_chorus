      module mod_structure

      use mod_common

      implicit none

!------------------------------------------------------------

      type structure_fils  
 
        type(structure_maille), pointer :: ptr => null()

      end type structure_fils

!------------------------------------------------------------

      type structure_support
 
        type(structure_maille), pointer :: suivant => null()
        type(structure_maille), pointer :: precedent => null()

      end type structure_support

!------------------------------------------------------------	

      type structure_maille 

! .....est-ce que la maille appartient a l'arbre ?
        logical :: arbre

! .....est-ce que la maille est reel (true) ou feuille fictive (false)
        logical :: feuille

! .....niveau de la maille
        integer :: maille_niv 

! .....valeur du seuil pour la maille
        real(kind=8) :: eps

! .....valeur du detail retenu pour la multiresolution
        real(kind=8) :: detail

! .....indice de la maille dans les directions (ndim)
        integer, dimension(ndim) :: maille_i 

! .....coordonnees et normale de la maille dans les directions (ndim)
        real(kind=8), dimension(ndim) :: x

! .....pas d'espace dans les directions (ndim)
        real(kind=8), dimension(ndim) :: dx

! .....valeur des details par variable (nvar) par rapport aux fils (2**ndim)
        real(kind=8), dimension(nvar, 2**ndim) :: det

! .....flux Euler a gauche et a droite par variable (nvar)
        real(kind=8), dimension(nvar) :: fmp_gch, fmp_dt

! .....flux visqueux Lax-Wendroff a gauche et a droite par variable (nvar)
!      et  par direction (ndim)
        real(kind=8), dimension(nvar, ndim) :: fvisc_gch, fvisc_dt

! .....variables conservatives au temps n*deltat
        real(kind=8), dimension(nvar) :: u_ndt

! .....variables conservatives au temps intermediaires pour M-C
        real(kind=8), dimension(nvar) :: u0

! .....variables conservatives au temps (n+1)*deltat
        real(kind=8), dimension(nvar) :: u
        real(kind=8), dimension(nvar) :: unp1

! .....pointer sur le pere de la maille
        type(structure_maille), pointer :: pere => null()

! .....pointer sur les fils (2**ndim) de la maille
        type(structure_fils), dimension(2**ndim) :: fils

! .....pointer sur les voisins par direction (ndim)
        type(structure_support), dimension(ndim) :: support

      end type structure_maille

!------------------------------------------------------------	

! .....pointer sur les mailles suivant indices structures (2**niv_max-1)**ndim
      TYPE(structure_fils), Pointer, DIMENSION(:) :: hash_table

! .....Liste des feuilles reelles ==> pointer sur la structure maille
      TYPE(structure_fils), Pointer, DIMENSION(:) :: feuille_reelle

!------------------------------------------------------------	

      end module mod_structure
