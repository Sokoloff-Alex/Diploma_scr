% test_3

PointA = [1.0 1.5];
PointB = [1.5 2.0];
PointC = [2.0 1.5];

VelA = [0.00  0.00];
VelB = [0.01  0.00];
VelC = [0.00  0.00];

clc
Strain = getStrain5(PointA, PointB, PointC, VelA, VelB, VelC)
