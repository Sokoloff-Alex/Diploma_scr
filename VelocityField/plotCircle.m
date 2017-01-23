function [pl] = plotCircle(lat, long, style)

Re = 6371*1000; %[km]
a = [1:360]';

r = Re*[cosd(a), sind(a), zeros(size(a))];

lat = lat-90;
long = long-90;

R1 = [1     0          0
      0  cosd(-long) sind(-long)
      0 -sind(-long) cosd(-long)];

R3 = [ cosd(-lat) sind(-lat)  0
      -sind(-lat) cosd(-lat)  0
          0         0       1];

r2 = zeros(length(a),3);
  
for i = 1:length(a)
    r2(i,:) = (R3*(R1*r(i,:)'))';
end
  
% figure(1)
% Earth_coast(3)
% plot3(r(:,1), r(:,2), r(:,3),'.-g')
plot3(r2(:,1),r2(:,2),r2(:,3),style)

end  
  
  