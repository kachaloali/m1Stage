#!bin/bash
#on renome tous les fichier .cif en .pdb
cd ./PDB/;
for ifile in *.$1;do mv $ifile ${ifile%$1}$2 ;done
cd ../
