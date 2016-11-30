% test_3

PointA = [1.0 1.5];
PointB = [1.5 2.0];
PointC = [2.0 1.5];

VelA = [0.00  0.00];
VelB = [0.00  0.01];
VelC = [0.00  0.00];

clc
format shortG
Strain = getStrain(PointA, PointB, PointC, VelA, VelB, VelC)
