#!/bin/bash
../PLM -TSfile Flat.q  -bilayerThickness 3.8 -rescalefactor 1.5
../PCG -dts point -str input.str -Bondlength 0.2 -LLIB Martini2.LIB
gmx editconf -f output.gro -o output.gro -box 22.50000  15.00000  15 -c
rm Flat.top
cat ./ff/include.txt output.top >> Flat.top
gmx grompp -f minimization.mdp -c output.gro -o EM -p Flat.top -maxwarn 2 -r output.gro
gmx mdrun -s EM -deffnm EM -v

gmx grompp -f EQ.mdp -c EM.gro -o EQ -p Flat.top -maxwarn 2 -r EM.gro
gmx mdrun -s EQ -deffnm EQ -v

../SOL -in EQ.gro -tem water.gro -o SOL.gro -ion 208 0
#gmx solvate -cp EQ.gro -cs water.gro -o SOL.gro  -radius 0.3 -p Flat.top
cat info.txt >> Flat.top

gmx grompp -f minimization.mdp -c SOL.gro -o EM -p Flat.top -maxwarn 2 -r SOL.gro
gmx mdrun -s EM -deffnm EM -v



gmx grompp -f EQS.mdp -c EM.gro -o EQ -p Flat.top -maxwarn 2 -r EM
gmx mdrun -s EQ -deffnm EQ -v

gmx grompp -f EQ.mdp -c EQ.gro -o EQ -p Flat.top -maxwarn 2 -r EQ

gmx genion -s EQ.tpr -p Flat.top  -neutral -o ION   -r EQ # -conc 0.05

gmx grompp -f EQS.mdp -c ION.gro -o EQ -p Flat.top -maxwarn 2 -r EQ
gmx mdrun -s EQ -deffnm EQ -v

gmx grompp -f PROD.mdp -c EQ.gro -o PROD -p Flat.top -maxwarn 2 -r EQ
gmx mdrun -s PROD -deffnm PROD -v



