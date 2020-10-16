
#!/bin/bash



../../Documents/TS2CG.1.1/TS2CG1.1/PLM -TSfile $1 -bilayerThickness 3.8 -rescalefactor 1.5

../../Documents/TS2CG.1.1/TS2CG1.1/PCG -dts point -str input.str -Bondlength 0.25 -WallBName WL -WallH 0.3  -Wall
