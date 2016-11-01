function [R_E2000, V_E2000] = ITRF2ETRF(R_I2008, V_I2008)
% Transformation of Velocity vector from  ITRF to ETRF  (substracting the plate motion)
% from ITRF2008 to ETRF2000
% checked at http://epncb.oma.be/_productsservices/coord_trans/index.php
% todo: 
% ITRF2008 -> ITRF2005 -> ITRF200 -> ETRF2000
% coordinates and velocities at epoch 2005.0
% parameters at epoch 2000.0
%
% Alexandr Sokolov, KEG
% 21.01.2016
    Tc = 2005.0; % central epoch

     R_I2000 = zeros(size(R_I2008));
     V_I2000 = zeros(size(V_I2008));

     R_E2000 = zeros(size(R_I2008));
     V_E2000 = zeros(size(V_I2008));
    
    mas2rad = 2*pi / (360 * 60 * 60 * 1000) ; % Miliarcsecond to radians % conv mas/y -> [rad/y
    
    % From table 5 # http://etrs89.ensg.ign.fr/memo-V8.pdf 
    % for epoch 2000.0  

     %%  ITRF2008 > ITRF2000 , ok !
%   transformation parameters from 
%   http://itrf.ensg.ign.fr/doc_ITRF/Transfo-ITRF2008_ITRFs.txt

    Rxyz    = [0 0 0];
    Rdotxyz = [0 0 0];
    
    D =    1.34 * 10^-9;
    Ddot = 0.08 * 10^-9;
    D = D + Ddot*(Tc - 2000);
    
    T =    [ -0.0019 -0.0017 -0.0105] ; % [m], ITRF2008 > ITRF2000 at 2000.0
    Tdot = [  0.0001  0.0001 -0.0018] ; % [m/yr]   
    T = T + Tdot*(Tc-2000); % propagate translation to the central epoch (here,  = 2005.0 by deafult)
    
    for i = 1:size(R_I2008,1)
        R_I2000(i,:) = R_I2008(i,:) + T    + D   * R_I2008(i,:);
        V_I2000(i,:) = V_I2008(i,:) + Tdot + Ddot* R_I2008(i,:);
    end 
   
%%   ITRF2000 > ETRF2000 ; ok !
    Rdot1 =  0.081 * mas2rad; 
    Rdot2 =  0.490 * mas2rad;
    Rdot3 = -0.792 * mas2rad;     
    Rdot = [  0    -Rdot3    Rdot2
            Rdot3     0     -Rdot1
           -Rdot2   Rdot1      0  ] ; % [rad/yr]   
    D    = 0;
    Ddot = 0;
   
    T =    [54.0   51.0  -48.0] * 0.001; % [m], ITRF2000 > ETRF2000
    Tdot = [ 0.0    0.0    0.0] * 0.001; % [m/yr]   
    T = T + Tdot*(Tc-2000); % propagate translation to the ecentral epoch (here,  = 2005.0 by deafult)

    for i = 1:size(R_I2000,1)
        R_E2000(i,:) = R_I2000(i,:) + T    + (Rdot * R_I2000(i,:)')'*(Tc-1989);
        V_E2000(i,:) = V_I2000(i,:) + Tdot + (Rdot * R_I2000(i,:)')';
    end
%     
% %%   ITRF2008 > ETRF2000 ; vel ok, crd not
% % 
%     R1 =  0.891 * mas2rad;
%     R2 =  5.390 * mas2rad;
%     R3 = -8.712 * mas2rad;
%     
%     R = [ 0   -R3    R2
%           R3   0    -R1
%          -R2   R1    0 ]; % [rad/yr]  
%     
%     Rdot1 =  0.081 * mas2rad; 
%     Rdot2 =  0.490 * mas2rad;
%     Rdot3 = -0.792 * mas2rad; 
%     Rdot = [  0    -Rdot3    Rdot2
%             Rdot3     0     -Rdot1
%            -Rdot2   Rdot1      0  ]; % [rad/yr]  
%        
% %     R = R + Rdot*(Tc-2000)
%        
%     D    = 1.34 * 10^-9;
%     Ddot = 0.08 * 10^-9;
%     D = D + Ddot*(Tc-2000);
%    
%     T =    [52.1  49.3  -58.5] * 0.001; % [m], ITRF2000 > ETRF2000
%     Tdot = [ 0.1   0.1   -1.8] * 0.001; % [m/yr]   
%     T = T + Tdot*(Tc-2000); % propagate translation to the ecentral epoch (here,  = 2005.0 by deafult)
% 
%     % equation for coordinates is not valid since the parameters are not
%     % simply the sum of two transformations
%     
%     for i = 1:size(R_I2008,1)
%         R_E2000(i,:) = R_I2008(i,:) + T +    D    * R_I2008(i,:) + (Rdot * R_I2008(i,:)')'*(Tc-1989);
%         V_E2000(i,:) = V_I2008(i,:) + Tdot + Ddot * R_I2008(i,:) + (Rdot * R_I2008(i,:)')';
%     end
% 
end