FC = ifort
LD = ifort

FFLAGS = -fpp -O3 -r8 -Wp,-DVISCCST,-DNbdim=2,-DOS7,-DCA,-DCav2
LDFLAGS = -fpp -r8 -Wp,-DVISCCST,-DNbdim=2,-DOS7,-DCA,-DCav2

## Tags :
## 1D => Nbdim=1, 2D => Nbdim=2 ou 3D => Nbdim=3
## TAC => TUBE À CHOC 1D,2D ou 3D
## ODC => Onde de CHOC 2D ou 3D
## CTI => Convection Tourbillon Isentropique 2D
## ODCS => Onde De Choc Stationnaire 1D,2D ou 3D
## SDT => Spot de Température, 1D ou 2D
## TDT => Tourbillon de Taylor 2D
## ICT1D => Interaction Choc-Turbulence 1D
## Cav1, Cav2 => Cavité entrainée 2D voir init_sol.f90

## OS7, TVDMP, ROE
## AC,CA,CAC,ACA

LINK = -c

LIBS =

.SUFFIXES: .f90 .o

.f90.o:
	$(FC) $(FFLAGS) -c $*.f90

OBJET = mod_common.o mod_structure.o mod_recherche.o \
	mod_initial.o mod_fonction.o mod_arbre.o mod_print.o mod_euler.o mod_save_vtk.o

mr_chorus.x: mr_chorus.o $(OBJET)
	$(LD) $(LDFLAGS) -o mr_chorus.x mr_chorus.o $(OBJET) $(LIBS)

mod_common.o : mod_common.f90
	$(FC) $(FFLAGS) $(LINK) mod_common.f90

mod_structure.o : mod_structure.f90 mod_common.o
	$(FC) $(FFLAGS) $(LINK) mod_structure.f90

mod_initial.o : mod_initial.f90 lect_don.f90 chainage.f90 support_racine.f90 \
init_sol.f90 lect_reprise.f90 recherche.f90 \
init_flags.f90 multimesh.f90 mod_common.o mod_structure.o mod_arbre.o
	$(FC) $(FFLAGS) $(LINK) mod_initial.f90

mod_fonction.o : mod_fonction.f90 \
liberation.f90 viscosite.f90 gradient.f90 mod_common.o mod_structure.o
	$(FC) $(FFLAGS) $(LINK) mod_fonction.f90

mod_recherche.o : mod_recherche.f90 recherche_support.f90 \
cal_indice.f90 creation.f90  mod_common.o mod_structure.o
	$(FC) $(FFLAGS) $(LINK) mod_recherche.f90

mod_arbre.o : mod_arbre.f90 encoding_moy.f90 encoding_det.f90 \
decoding_moy.f90 seuillage.f90 graduation.f90 graduation_local.f90 \
elagage.f90 chainage_support.f90 liste_feuille.f90 feuille_fictive.f90 \
feuille_fictive_visc.f90 valeur_feuille_fictive.f90 \
mod_common.o mod_structure.o mod_fonction.o mod_recherche.o
	$(FC) $(FFLAGS) $(LINK) mod_arbre.f90

mod_print.o : mod_print.f90 print_leaves.f90 print_reprise.f90 \
mod_common.o mod_structure.o mod_fonction.o
	$(FC) $(FFLAGS) $(LINK) mod_print.f90

mod_euler.o : mod_euler.f90 pas_de_temps_sp.f90 \
undt_update.f90 u0_update.f90 u_update_acc.f90 u_update_conv.f90\
flux_euler.f90 flux_visc.f90 \
flux_osmp7_conv.f90 flux_osmp7_acoustic.f90 \
vcalir.f90 vcalr.f90 vcalir_acc.f90 osmp7_conv.f90 osmp7_acoustic.f90 \
integration_euler.f90  integration_euler_acc.f90 \
integration1_visqueux.f90 integration2_visqueux.f90 \
mise_a_jour.f90 caldeb_flux_euler.f90 calfin_flux_euler.f90 \
caldeb_flux_euler_ac.f90 calfin_flux_euler_ac.f90 \
caldeb_flux_visc.f90 calfin_flux_visc.f90 residu.f90 \
mod_common.o mod_structure.o mod_recherche.o mod_fonction.o
	$(FC) $(FFLAGS) $(LINK) mod_euler.f90

mod_save_vtk.o : mod_save_vtk.f90  \
mod_common.o mod_structure.o
	$(FC) $(FFLAGS) $(LINK) mod_save_vtk.f90

mr_chorus.o : mr_chorus.f90 mod_common.o mod_structure.o \
mod_initial.o mod_arbre.o mod_recherche.o \
mod_print.o mod_fonction.o mod_euler.o \
mod_save_vtk.o
	$(FC) $(FFLAGS) $(LINK) mr_chorus.f90

clean:
	rm -f *.o *.a *.x *.mod

cdat:
	rm I*.dat T*.dat

again: clean all

##depend:
##	f90 *.f90 -MDupdate .depends
##
##.depends:
##	touch .depends
##
##include .depends

