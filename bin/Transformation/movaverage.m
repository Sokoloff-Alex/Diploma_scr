function [vector_ave] = movaverage(vector, win_size)
% computes centered moving average

if size(vector,1) == 1 && size(vector,2) ~= 1
    vector = vector';
end

vector_ave = zeros(size(vector));
shift = (win_size-1)/2;
vector2 = [nan(shift,1); vector; nan(shift,1)];

for i = 1:length(vector)   
    vector_ave(i) = nanmean(vector2(i:(i+win_size-1)));
end

if size(vector,1) ~= 1 && size(vector,2) == 1
    vector_ave = vector_ave';
end

end