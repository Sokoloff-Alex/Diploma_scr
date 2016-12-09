function [LongGrid, LatGrid, V_def3, rmsFit, V_SigPred ] = run_Collocation(long, lat, Venu, CovVenu, LongLim, LatLim, step, Max_Dist, nLim )
% run Collocation

% LSC, % allocation domain for grid points
tic
range = 1:length(lat);


nLon = ceil(( LongLim(2) - LongLim(1) )/step) +1;
nLat = ceil(( LatLim(2)  - LatLim(1)  ) / step) +1;
n = nLat * nLon ;

% Dummies
V_def3    = zeros(n,3);
rmsFit    = zeros(n,6);
V_SigPred = zeros(n,3);
LatGrid   = zeros(n,1);
LongGrid  = zeros(n,1);

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
%             Venu = [Ve_res(sel), Vn_res(sel), Vu_res(sel)]*1000;
        [V_pred, rmsFitting, V_noise_pred] = solve_WLSC3(iLat, iLong, lat(sel), long(sel), Venu(sel,:)*1000 ,CovVenuSel,'exp1', '-v', 'bias', 'tail 0', 'no corr', 'filter');   % WLSC 
        V_def3(p,:) = V_pred/1000;
        rmsFit(p,:)   = rmsFitting/1000; % [mm/yr]
        V_SigPred(p,:) = V_noise_pred;
        LatGrid(p,1)  = iLat;
        LongGrid(p,1) = iLong;
    end
end

toc

end