function [V_vertical, V_horizonal] = VerticalProjection(X, Y, Z, Vx, Vy, Vz)
% compute vertical projection of velocity vector

len = length(X);
V_vertical = zeros(len,3);
V_horizonal = zeros(len,3);

for i = 1:len
    
    e_radial(i,:) = [X(i),Y(i),Z(i)]./(norm([X(i),Y(i),Z(i)]));
    e_v(i,:) = [Vx(i), Vy(i), Vz(i)]./(norm([Vx(i), Vy(i), Vz(i)]));
    
    e_o(i,:) = cross(e_radial(i,:), e_v(i,:));
    e_h(i,:) = cross(e_o(i,:), e_radial(i,:));
        
    theta = acosd(e_radial(i,:) * e_v(i,:)');
    
    V_vertical(i,:) = e_radial(i,:).* norm([Vx(i), Vy(i), Vz(i)]) * cosd(theta);
    V_horizonal(i,:) = e_h(i,:)    .* norm([Vx(i), Vy(i), Vz(i)]) * sind(theta);
    
end

end