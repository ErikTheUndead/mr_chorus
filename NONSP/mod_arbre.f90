
	module mod_arbre

	use mod_common
	use mod_structure
	use mod_fonction
	use mod_recherche

	contains

	include "encoding_moy.f90"
	include "encoding_det.f90"
	include "seuillage.f90"
	include "graduation.f90"
	include "graduation_local.f90"
	include "elagage.f90"
	include "chainage_support.f90"
	include "liste_feuille.f90"
	include "feuille_fictive.f90"
	include "feuille_fictive_visc.f90"
        include "valeur_feuille_fictive.f90"
	include "decoding_moy.f90"

	end module mod_arbre
