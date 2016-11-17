% test OrthodromeArcLength

clc

%%
clc
p = length(lat);
arc = zeros(p,p);

for i = 1:p 
    stack = OrthodromeArcLength(lat(i), long(i), lat, long);   
    for j = 1:p
       if ~isreal(stack(j))
          disp(['Error :: SITES:: i = ',num2str(i),...
                        ' j=',num2str(j),...    
                        ' lat1=', num2str(lat(i)), ...
                        ' long1=',num2str(long(i)), ...
                        ' lat2=', num2str(lat(i)), ...
                        ' long1=',num2str(long(i))]); 
       end  
    end
    arc(:,i) = stack;
    if arc(i,i) ~= 0 
        disp(['Error, must be 0 Arc !!!, i = ',num2str(i), ', arc=',num2str(arc(i,i)) ])
    end 
end

%% compare with Haversine equation
clc
for i = 1:p
   for j = 1:p
       arc1 = OrthodromeArcLength(lat(i), long(i), lat(j), long(j));
       arc2 = greatcircleArc(lat(i), long(i), lat(j), long(j));
       Error = arc1 - arc2;
       if Error ~= 0
          disp(['Error: arc1 ~= arc2, Erorr = ', num2str(Error)]) 
       end
   end
end
















