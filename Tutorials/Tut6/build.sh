#!/bin/csh


source /usr/local/gromacs-2020/bin/GMXRC
set ts2cg = '/coarse/weria/Projects/TS_Martini_Development/TS2CG1.1'



$ts2cg/PCG  -str $1 -Bondlength 0.3 -LLIB  Martini2.LIB  -function analytical_shape


