% test_kriging

close all
clc

Outliers = {'HELM','WIEN','FERR','FERH','OGAG', ...
            'ROHR','BIWI','BI2I','MANS','FFMJ','MOGN','WLBH', ...
            'TRF2','KRBG','OBE4','WT21','HKBL','PATK','PAT2', ...
            'VAUC','HRIE','ENTZ','KTZ2','BZRG'};
     
iiOut = selectRange(names, Outliers);
ii = 1:length(Vu_res);
iiSel = setdiff(ii, iiOut);

x = wrapTo180(long(iiSel));
y = lat(iiSel);
z = Vu_res(iiSel)*1000;
%% Combine Predicted and Original
% 
x = [LongGrid; wrapTo180(long(iiSel))];
y = [LatGrid; lat(iiSel)];
z = [V_def3(:,3); Vu_res(iiSel)] *1000;


%% WLCS only
x = LongGrid;
y = LatGrid;
z = V_def3(:,3) *1000;

%%
% interpolate default
[X, Y] = meshgrid(0:0.05:17, 42:0.05:50);
y = y * (4/3);
Y = Y * (4/3);
vq = griddata(x, y, z, X, Y, 'natural');
clc

try
    close (fig1)
end
fig1 = figure(1);
subplot(2,2,1)
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
contour(X,Y,Zhat,[-2:0.5:3],'LineWidth',2)
% axis image
% quiver(LongGrid, LatGrid, zeros(size(V_def3(:,3))), V_def3(:,3)*100, 0, 'b')
% quiver(x, y, zeros(size(z)), z/1000*200, 0, 'b')
quiver(long,     lat,     zeros(size(Vu_res)),      Vu_res*100,      0, 'k', 'LineWidth',1)

% quiver3(x, y, abs(z), zeros(size(z)),zeros(size(z)), z,'b')
% text(x(z >= 0),y(z >= 0),2*z(z >= 0),names(iiSel(z >= 0)))
% text(x(z < 0),y(z < 0),abs(z(z<0)),names(iiSel(z < 0)))
text(long(iiOut), lat(iiOut), names(iiOut), 'Color', 'r')
text(long(iiSel), lat(iiSel), names(iiSel), 'Color', 'b')
text(long, lat, num2str(Vu_res*1000, '%4.1f'), 'HorizontalAlignment', 'right');
% text(LongGrid, LatGrid, num2str(V_def3(:,3)*1000, '%4.1f'), 'HorizontalAlignment', 'right');
colorbar
title('kriging predictions')
xlim([2 17])
ylim([43 50])
zlim([-10 10])

%% Save for GMT

a = grid2vector(X);
b = grid2vector(Y);
c = grid2vector(Zhat);

writeVelocityFieldVertical2GMT(a,b,c, '~/Alpen_Check/MAP/VelocityField/Vel_up_Kriging2.txt');
% writeVelocityFieldVertical2GMT(long(iiSel),lat(iiSel),Vu_res(iiSel)*1000, '~/Alpen_Check/MAP/VelocityField/Vu_res.txt');





