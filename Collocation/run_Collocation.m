function [LongGrid_stack, LatGrid_stack, V_pred_stack, rmsFit_stack, V_SigPred_stack ] = run_Collocation(long, lat, Venu, CovVenu, LongLim, LatLim, step, Max_Dist, nLim, varargin)
% run Collocation

% LSC, % allocation domain for grid points
tic;

flags = varargin(:);

range = 1:length(lat);


nLon = ceil(( LongLim(2) - LongLim(1) ) / step) +1;
nLat = ceil(( LatLim(2)  - LatLim(1)  ) / step) +1;
n = nLat * nLon ;

% Dummies
V_pred_stak     = zeros(n,3);
rmsFit_stack    = zeros(n,6);
V_SigPred_stack = zeros(n,3);
LatGrid_stack   = zeros(n,1);
LongGrid_stack  = zeros(n,1);

p = 0;
for iLong = LongLim(1):step:LongLim(2)
    for iLat = LatLim(1):step:LatLim(2)
        arc = greatcircleArc(iLat, iLong, lat, long) * 111 ; % km
        sel = range(arc < Max_Dist);    
        if length(sel) < nLim % add more stations
            add = sort(arc);
            add = add(1:nLim); % ad 2 sites, ast is prop already included
            iadd = ismember(arc,add);
            inew = range(iadd);
            sel = unique(sort([sel, inew]));
        end
        p = p + 1; % point Number
        CovVenuSel = extractCovariance(CovVenu, sel, [1 2 3], 'split');
        % WLSC
        [V_pred, rmsFit, V_noise_pred] = solve_WLSC3(iLat, iLong, lat(sel), long(sel), Venu(sel,:)*1000 ,CovVenuSel, flags);
        V_pred_stack(p,:)    = V_pred/1000;
        rmsFit_stack(p,:)    = rmsFit/(1000); % [mm/yr]
        V_SigPred_stack(p,:) = V_noise_pred;
        LatGrid_stack(p,1)   = iLat;
        LongGrid_stack(p,1)  = iLong;
    end
end

toc

end