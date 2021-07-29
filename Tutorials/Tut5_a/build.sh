#!/bin/csh


source /usr/local/gromacs-2020/bin/GMXRC
set ts2cg = '/coarse/weria/Projects/TS_Martini_Development/TS2CG1.1'




$ts2cg/PLM -TSfile   ./$1   -bilayerThickness 4 -rescalefactor 6 6 6 
$ts2cg/PCG   -str input.str -Bondlength 0.2 -LLIB Martini2.LIB 



