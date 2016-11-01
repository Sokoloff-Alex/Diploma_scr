function [V_res_blocks, V_blocks, Omega_Blocks] = remove_block_motion(CRD, VEL, Sets)
% 
% function to remove motion in each block individually
% input : CRD, VEL(plate motion removed)
%       : cellArray of Sets for each block
%
% Alexandr Sokolov, KEG
% 31.10.2016
    
    Omega_Eur = [55.9533, -97.4134,   2.6364e-07 ]';
    V_res_blocks = VEL;
    Omega_Blocks = zeros(length(Sets),3);

    for i = 1:length(Sets)
        iSet = Sets{i};
        % estimate Euler pole for a block
        [Omega_Blocks(i,:)] = plate_motion(Omega_Eur, CRD(iSet,:),VEL(iSet,:), 500);

        % remove block motion
        [V_res_blocks(iSet,:), V_blocks(iSet,:)] = remove_plate_motion(CRD(iSet,:), VEL(iSet,:), Omega_Blocks(i,:));
        
        % check for 3 sigma
        [Ve_res_bl, Vn_res_bl, Vu_res_bl] = XYZ2ENU(CRD,V_res_blocks); 
        disp(['iSet # ', num2str(i), ' ; ::: # of sites:', num2str(length(iSet)) ])
        disp(['std Ve_res_bl = ', num2str(std(Ve_res_bl)*1000), '[mm/yr]' ])
        disp(['std Vn_res_bl = ', num2str(std(Vn_res_bl)*1000), '[mm/yr]' ])
        
    end
end