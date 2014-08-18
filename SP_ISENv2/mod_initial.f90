      module mod_initial

      use mod_common
      use mod_structure
      use mod_arbre


      contains

!      include "lect_don.f90"
#include <lect_don.f90>
      include "chainage.f90"
      include "support_racine.f90"
      
#include <init_sol.f90>
!      include "init_sol.f90"

      include "lect_reprise.f90"
      include "recherche.f90"

      include "init_flags.f90"
      include "multimesh.f90"

      end module mod_initial
