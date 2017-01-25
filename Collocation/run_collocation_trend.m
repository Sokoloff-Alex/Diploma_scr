function [V_trend_stack, V_sigma_trend_stack] = run_collocation_trend(long, lat, V_enu_res, CovVenu, Cov_scale, Max_Dist, nLim, varargin)
% function ro run collocation to estimate the trend model.

flags = varargin(:)

V_trend_stack = zeros(size(V_enu_res));
V_sigma_trend_stack = zeros(size(V_enu_res));


tic
for i = 1:length(long) 
    [LongGrid, LatGrid, V_trend, rmsFit, V_SigTrend] = run_Collocation(long, lat, V_enu_res, CovVenu, Cov_scale, [long(i) long(i)], [lat(i) lat(i)], 1, Max_Dist, nLim, 'exp1', '-v', 'bias', 'tail 0', 'no corr', 'filter');
    V_trend_stack(i,:) = V_trend;
    V_sigma_trend_stack(i,:) = V_SigTrend;
end

toc

end