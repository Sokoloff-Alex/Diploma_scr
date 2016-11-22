function [] = ellipce_2D(sigmas, azimuth, origin, scale, varargin)
% ellipses 2D 

if isempty(varargin)
    clr = [0.5 0.5 0.5];
else
    clr = varargin{1};
end
    
    
if sigmas(1) == 0
   sigmas(1) = 0.00001; 
end

if sigmas(2) == 0
   sigmas(2) = 0.00001; 
end

%%
NP = 50;
xy = zeros(NP, 2);
alpha  = 2*pi/NP*(0:NP);
circle = [cos(alpha);sin(alpha)]; 
P = [sigmas(1)      0; 
        0      sigmas(2)]*scale;
ellip = P*circle;
X = ellip(1,:);
Y = ellip(2,:);

X = [X, X(1)];
Y = [Y, Y(1)];

% Rotate ellipse
% Orientation of the error ellipse: azimuth of the principal axis
% counted positive in East (in degrees)
if P(1,1) < P(2,2)
   azimuth = azimuth+90; 
end

R = [cosd(azimuth) -sind(azimuth)
     sind(azimuth)  cosd(azimuth)];

 for i=1:NP+1
     xy(i,:) = R*[X(i); Y(i)] + origin';
 end

plot(xy(:,1), xy(:,2),'-', 'Color',clr);
hold on

end