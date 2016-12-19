function [LongK_Stack, LatK_Stack, VuK_Stack , fig1, fig2] = runKriging(LongGrid, LatGrid, V_def, long ,lat,  Vu_res, iiSel) 
% function to run Kriging 


% WLCS only
x = LongGrid;
y = LatGrid;
z = V_def(:,3) *1000;

%% interpolate by default Matlab interpolator
[X, Y] = meshgrid(min(LongGrid):0.1:max(LongGrid), min(LatGrid):0.1:max(LatGrid));
y = y * (4/3);
Y = Y * (4/3);
vq = griddata(x, y, z, X, Y, 'natural');
clc

%% do Kriging
try
    close (fig1)
end
fig1 = figure(1);
subplot(2,2,1)
colormap('jet')
% scatter3(x,y,z); axis image; axis xy
hold on
Earth_coast(2)
% plot(x,y,'.k')
% quiver(x, y, zeros(size(z)), z,'b')

pcolor(X,Y *(3/4),vq);
% text(x,y,names )
shading interp
xlim([4 14])
ylim([44 48])
colorbar
title('random field with sampling locations')


% calculate the sample variogram
v = variogram([x y],z,'plotit',false,'maxdist',5);
% and fit a spherical variogram
subplot(2,2,2)
[dum,dum,dum,vstruct] = variogramfit(v.distance,v.val,[],[],[],'model','stable');
title('variogram')

% now use the sampled locations in a kriging
[Zhat,Zvar] = kriging(vstruct,x,y,z,X,Y,50);
subplot(2,2,3)
hold on
axis xy
y = y * (3/4);
Y = Y * (3/4);
pcolor(X(1,:),Y(:,1),Zhat); %
shading interp
% quiver(x, y, zeros(size(z)), z,'b')
contour(X,Y,Zhat)
% axis image

colorbar
title('kriging predictions')
subplot(2,2,4)
hold on
% contour(X,Y,Zvar); %axis image
contour(X,Y,Zhat)
axis xy
% title('kriging variance')#
% title('kriging prediction'

%%
try
    close(fig2)
    clc
end
fig2 = figure(2);
hold on
grid on
% etopo_fig = showETOPO(ETOPO_Alps.Etopo_Europe, ETOPO_Alps.refvec_Etopo);
Earth_coast(2)
h = pcolor(X(1,:),Y(:,1),Zhat); %
shading interp
set(h,'facealpha',.5)
colormap('jet')
contour(X,Y,Zhat,[-2:0.5:-.5, 0.5:0.5:3],'LineWidth',2)
% contour(X,Y,Zhat,[-2:0.5:3],'LineWidth',2)

% axis image
% quiver(LongGrid, LatGrid, zeros(size(V_def3(:,3))), V_def3(:,3)*100, 0, 'b')
% quiver(x, y, zeros(size(z)), z/1000*200, 0, 'b')
quiver(long,     lat,     zeros(size(Vu_res)),      Vu_res*100,      0, 'k', 'LineWidth',1)

% quiver3(x, y, abs(z), zeros(size(z)),zeros(size(z)), z,'b')
% text(x(z >= 0),y(z >= 0),2*z(z >= 0),names(iiSel(z >= 0)))
% text(x(z < 0),y(z < 0),abs(z(z<0)),names(iiSel(z < 0)))
% text(long(iiOut), lat(iiOut), names(iiOut), 'Color', 'r')
% text(long(iiSel), lat(iiSel), names(iiSel), 'Color', 'b')
plot(long(iiSel), lat(iiSel), '.')
text(long, lat, num2str(Vu_res*1000, '%4.1f'), 'HorizontalAlignment', 'right');
% text(LongGrid, LatGrid, num2str(V_def(:,3)*1000, '%4.1f'), 'HorizontalAlignment', 'right');
colorbar
title('kriging predictions')
xlim([min(LongGrid) max(LongGrid)])
ylim([min(LatGrid)   max(LatGrid)])
% zlim([-10 10]) % [mm/yr]

%% Save for GMT

LongK_Stack = grid2vector(X);
LatK_Stack = grid2vector(Y);
VuK_Stack = grid2vector(Zhat);

% writeVelocityFieldVertical2GMT(LongK_Stack,LatK_Stack,VuK_Stack, '~/Alpen_Check/MAP/VelocityField/Vel_up_Kriging2.txt');
% writeVelocityFieldVertical2GMT(long(iiSel),lat(iiSel),Vu_res(iiSel)*1000, '~/Alpen_Check/MAP/VelocityField/Vu_res.txt');


end


