function [L1, L2] = netRotation(x1, y1, x2, y2)
% estimate net rotaton
%
% Alexandr Sokolov , DGFI
% 07.11.2016

A = [1 x1 y1  0  0 0;
        0  0  0 1 x1 y1];
    
X = (A'*A)^-1*A'*[x2; y2];

a0 =  X(1);
a1 =  X(2);
a2 =  X(3);
b0 =  X(4);
b1 =  X(5);
b2 =  X(6);


% NetRotation
omega = atan((a2-b1)/(a1+b1));

%
r = a1*cos(omega) - b1*sin(omega);
s = a1*sin(omega) + b1*cos(omega);
t = a2*sin(omega) + b2*cos(omega);

% Strains
L1 = (r+t)/2 + 1/2*sqrt((r-t)^2 + 4*s) - 1;
L2 = (r+t)/2  -  1/2*sqrt((r-t)^2 + 4*s) - 1;

end