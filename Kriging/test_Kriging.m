% test_kriging

close all
clc

Out   = {'HELM', 'WIEN', 'FERR', 'FERH', 'OGAG', 'SOND', 'OBE2', ...
                'ROHR','BIWI','BI2I','MANS', 'FFMJ', 'MOGN', 'WLBH', 'TRF2'};
iiOut    = selectRange(names, Out);

ii = 1:length(Vu_res);

iiSel = setdiff(ii, iiOut);

%
[X, Y] = meshgrid(-4:0.25:18, 42:0.25:52);

x = wrapTo180(long(iiSel));
y = lat(iiSel);
z = Vu_res(iiSel)*1000;

%
close all
figure(1)
subplot(2,2,1)
% scatter3(x,y,z); axis image; axis xy
hold on
% plot(x,y,'.k')
quiver(x, y, zeros(size(z)), z,'b')
Flags_All_good = Selected;
vq = griddata(x, y, z, X, Y, 'natural');
pcolor(X,Y,vq);
% text(x,y,names )
shading interp
colorbar
title('random field with sampling locations')

% calculate the sample variogram
v = variogram([x y],z,'plotit',false,'maxdist',100);
% and fit a spherical variogram
subplot(2,2,2)
[dum,dum,dum,vstruct] = variogramfit(v.distance,v.val,[],[],[],'model','stable');
title('variogram')

% now use the sampled locations in a kriging
[Zhat,Zvar] = kriging(vstruct,x,y,z,X,Y);
subplot(2,2,3)
hold on
pcolor(X(1,:),Y(:,1),Zhat); %
shading interp
contour(X,Y,Zhat)
% axis image
axis xy
colorbar
title('kriging predictions')
subplot(2,2,4)
hold on

contour(X,Y,Zvar); %axis image
axis xy
colorbar
title('kriging variance')

