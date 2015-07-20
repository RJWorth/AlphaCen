#!/bin/bash
##############################################################################
# This script cleanly reruns the most recent Merc95

rm Out/*
rm Aei/*

gfortran -w -o merc Code/mercury6_2.f95 Code/drift.f95 Code/orbel.f95 Code/mal.f95 Code/mce.f95 Code/mco.f95 Code/mdt.f95 Code/mio.f95 Code/mfo.f95 Code/mxx.f95 Code/both.f95
./merc
gfortran -w -o elem Code/element6.f95 Code/e_sub.f95 Code/orbel.f95 Code/both.f95
./elem

#mv *.out $1
#mv *.dmp $1
mv *.tmp Out
mv *.aei Aei

