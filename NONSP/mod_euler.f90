      module mod_euler

      use mod_common
      use mod_structure
      use mod_recherche
      use mod_fonction

      contains

      include "pas_de_temps.f90"
      include "undt_update.f90"
      include "u0_update.f90"
      include "flux_euler.f90"
      include "flux_osmp7.f90"
      include "flux_visc.f90"

      include "vcalir.f90"
      include "vcalr.f90"
#include "osmp7.f90"

      include "integration_euler.f90"
      include "integration1_visqueux.f90"
      include "integration2_visqueux.f90"

      include "mise_a_jour.f90"

      include "caldeb_flux_euler.f90"
      include "calfin_flux_euler.f90"
      include "caldeb_flux_visc.f90"
      include "calfin_flux_visc.f90"

      include "residu.f90"

      end module mod_euler
