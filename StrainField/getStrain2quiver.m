function [Stack_dilat, Stack_compr] =  getStrain2quiver(Strain, sc)
% function for prepare tables (Stack) to plot by quiver
%
% input : Strain - stack of strain field parameters
%         sc     -  scale
%
% output: Stack_dilat, Stack_compr - sepatate tables for dilatiationa and
% compression, for show Normal Strain field
%
% command:
%   [Stack_dilat, Stack_compr] =  getStrain2quiver(Strain)
%
% Alexandr Sokolov, KEG
% 24.11.2016

Long  = Strain(:,1); 
Lat   = Strain(:,2);
L1    = Strain(:,3);
L2    = Strain(:,4);
Omega = Strain(:,5);

N = length(L1);

%% prepare qstack for quiver

nPosL1 = length(L1(L1 > 0));
nPosL2 = length(L2(L2 > 0));

Stack_L1_pos_a = NaN(nPosL1,4);
Stack_L1_pos_b = NaN(nPosL1,4);

Stack_L1_neg_a = NaN(N-nPosL1,4);
Stack_L1_neg_b = NaN(N-nPosL1,4);

Stack_L2_pos_a = NaN(nPosL2,4);
Stack_L2_pos_b = NaN(nPosL2,4);

Stack_L2_neg_a = NaN(N-nPosL2,4);
Stack_L2_neg_b = NaN(N-nPosL2,4);

counter_pos_L1 = 0;
counter_neg_L1 = 0;
counter_pos_L2 = 0;
counter_neg_L2 = 0;

for i = 1:size(Strain,1)
    % first axis
    dx_L1 = L1(i)*cosd(Omega(i));
    dy_L1 = L1(i)*sind(Omega(i));

    % second axis
    dx_L2 = L2(i)*cosd(Omega(i)+90);
    dy_L2 = L2(i)*sind(Omega(i)+90);

%   do L1
    if L1(i) > 0
        counter_pos_L1 = counter_pos_L1 + 1;
        Stack_L1_pos_a(counter_pos_L1,:) = [Long(i), Lat(i),  dx_L1*sc,  dy_L1*sc ];
        Stack_L1_pos_b(counter_pos_L1,:) = [Long(i), Lat(i), -dx_L1*sc, -dy_L1*sc ];      
    else 
        counter_neg_L1 = counter_neg_L1 + 1;
        Stack_L1_neg_a(counter_neg_L1,:) = [Long(i) + dx_L1*sc, Lat(i) + dy_L1*sc, -dx_L1*sc, -dy_L1*sc];
        Stack_L1_neg_b(counter_neg_L1,:) = [Long(i) - dx_L1*sc, Lat(i) - dy_L1*sc, +dx_L1*sc, +dy_L1*sc];
    end 
    
    %   do for L2
    if L2(i) > 0
        counter_pos_L2 = counter_pos_L2 + 1;
        Stack_L2_pos_a(counter_pos_L2,:) = [Long(i), Lat(i),  dx_L2*sc,  dy_L2*sc ];
        Stack_L2_pos_b(counter_pos_L2,:) = [Long(i), Lat(i), -dx_L2*sc, -dy_L2*sc ];      
    else 
        counter_neg_L2 = counter_neg_L2 + 1;
        Stack_L2_neg_a(counter_neg_L2,:) = [Long(i) + dx_L2*sc, Lat(i) + dy_L2*sc, -dx_L2*sc, -dy_L2*sc];
        Stack_L2_neg_b(counter_neg_L2,:) = [Long(i) - dx_L2*sc, Lat(i) - dy_L2*sc, +dx_L2*sc, +dy_L2*sc];
    end       
end

Stack_dilat = [Stack_L1_pos_a; Stack_L1_pos_b; Stack_L2_pos_a; Stack_L2_pos_b];
Stack_compr = [Stack_L1_neg_a; Stack_L1_neg_b; Stack_L2_neg_a; Stack_L2_neg_b];

end