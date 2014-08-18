#!/bin/bash
#PBS -N Choc_Spot
#PBS -l cput=72:00:00,mem=2600mb

set -x

DIR_PROG=$people/duhem/stage_2014_Limsi/MR_CHORUS_SPLITTING
DIR_DONNEES=$WORKDIR/
DIR_OUT=$WORKDIR/
#
mkdir -p $DIR_OUT

#--------------------------------------------
# On se place dans le TMPDIR de stokes
#--------------------------------------------
cd $TMPDIR

fget  $DIR_PROG/don_mrchorus.dat don_mrchorus.dat
fget  $DIR_PROG/mr_chorus.x mr_chorus.x
# fget  $DIR_DONNEES/save_mrchorus.dat save_mrchorus.dat

ls -lrt

#--------------------------------------------
echo "Execution sur stokes"
#--------------------------------------------
date
time ./mr_chorus.x

#--------------------------------------------
echo "Sauvegarde du fichier de resultats sur stokes"
#--------------------------------------------
ls -lrt

date

fput res_mrchorus.dat $DIR_OUT/res_mrchorus.dat

fput save_mrchorus.dat $DIR_OUT/save_mrchorus.dat

for i in $(ls Tree*.dat)
do
fput  $i $DIR_OUT/$i
done

for i in $(ls Int_OS7DIF2_*.dat)
do
fput  $i $DIR_OUT/$i
done

