% validation_swiss
%
% compare ALP_NET velocities with swisstopo results
%
%  In ETRF2000:
%  http://pnac.swisstopo.admin.ch/divers/dens_vel/ch08_hor.vel
%  http://pnac.swisstopo.admin.ch/divers/dens_vel/ch08_ver.vel
% 
%  in IGb08:
%  http://pnac.swisstopo.admin.ch/divers/dens_vel/ch08_hor_i08.vel
%  http://pnac.swisstopo.admin.ch/divers/dens_vel/ch08_ver_i08.vel
%
% ch08veri08(:,1) = table(table2array(ch08veri08(:,1)) - 360);
% ch08hori08(:,1) = table(table2array(ch08hori08(:,1)) - 360);
% ch08ver(:,1) = table(table2array(ch08ver(:,1)) - 360);
% ch08hor(:,1) = table(table2array(ch08hor(:,1)) - 360);

% ch08veri08.Properties.RowNames = table2array( ch08veri08(:,8) );
% ch08hori08.Properties.RowNames = table2array( ch08hori08(:,8) );
% ch08ver.Properties.RowNames = table2array( ch08ver(:,8) );
% ch08hor.Properties.RowNames = table2array( ch08hor(:,8) );

%%
% save('../../dat/ch08veri08.mat', 'ch08veri08')
% save('../../dat/ch08hori08.mat', 'ch08hori08')
% save('../../dat/ch08ver.mat', 'ch08ver')
% save('../../dat/ch08hor.mat', 'ch08hor')

%% load
load('../../dat/ch08veri08.mat')
load('../../dat/ch08hori08.mat')
load('../../dat/ch08ver.mat')
load('../../dat/ch08hor.mat')

%% remove plate motion


[V_res_xyz_false, V_pl_1] = remove_plate_motion(CRD(iiAlSt,:), VEL(iiAlSt,:), Omega_Eur);
[Ve_pl1, Vn_pl1, Vu_pl1] = XYZ2ENU(CRD(iiAlSt,:),V_pl_1);

swVu = table2array( ch08veri08(names_common,4) ) ;
swVe_res = table2array( ch08hori08(names_common,3) )- Ve_pl1*1000;
swVn_res = table2array( ch08hori08(names_common,4) ) - Vn_pl1*1000;


table_sw_res = table(ts.VarName1, ts.lat, swVe_res, swVn_res ,swVu, names_common, 'RowNames', names_common)


%%
names_alp = names;
alp_table = table(lat, long, Ve*1000, Vn*1000, Vu*1000, Ve_res*1000, Vn_res*1000, Vu_res*1000, names_alp, 'RowNames', names_alp)

%% compare velocities for common stations

names_swiss = ch08veri08.Properties.RowNames
names_common = intersect(names_swiss, names_alp)
r = 1:198;
iiAlSt = r(ismember(names, names_common))'

%% plot

% Defaults for this blog post
width =  7;     % Width in inches
height = 3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 9;      % Fontsize
lw =  1;      % LineWidth
msz = 3;       % MarkerSize

close all
% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);

close all
fig1 = figure(1);

subplot(1,3,1)
hold on;  axis equal
plot([16.8 22],[16.8 22],'-b')
Ve_ch_comm = ch08hori08{names_common,3};
Ve_al_comm = alp_table{names_common,3};
plot( Ve_ch_comm , Ve_al_comm , 'or','MarkerFaceColor','r') 
title(['V_E_a_s_t, mm/yr'])
% text(17, 21.5, [ 'mean(dV_e) = ', num2str(mean(Ve_ch_comm - Ve_al_comm), '%5.3f'), ' mm/yr'])
text(17, 21.5,   [ 'rms = ',  num2str(rms( Ve_ch_comm - Ve_al_comm), '%5.2f'), ' mm/yr'])
xlabel('Swisstopo')
ylabel('ALPNET')
xlim([16.8 22])
ylim([16.8 22])

subplot(1,3,2)
hold on; axis equal
plot([15 18],[15 18],'-b')
Vn_ch_comm = ch08hori08{names_common,4};
Vn_al_comm = alp_table{names_common,4};
plot( Vn_ch_comm , Vn_al_comm , 'or','MarkerFaceColor','r') 
title(['V_N_o_r_t_h, mm/yr'])
% text(15.2, 17.8, [ 'mean(dV_e) = ', num2str(mean(Vn_ch_comm - Vn_al_comm), '%5.3f'), ' mm/yr'])
text(15.2, 17.8,   [ 'rms = ',  num2str(rms( Vn_ch_comm - Vn_al_comm), '%5.2f'), ' mm/yr'])
xlabel('Swisstopo')
ylabel('ALPNET')
xlim([15 18])
ylim([15 18])

subplot(1,3,3)
hold on;  axis equal
plot([-1.8 2.5],[-1.8 2.5],'-b')
Vu_ch_comm = ch08veri08{names_common,4};
Vu_al_comm = alp_table{names_common,5};
plot( Vu_ch_comm , Vu_al_comm , 'or', 'MarkerFaceColor','r') 
title(['V_U_p, mm/yr'])
% text(-1.6, 2.2, [ 'mean(dV_u) = ', num2str(mean(Vu_ch_comm - Vu_al_comm), '%5.3f'), ' mm/yr'])
text(-1.6, 2.1, [ 'rms = ',  num2str(rms( Vu_ch_comm - Vu_al_comm), '%5.2f'), ' mm/yr'])
xlabel('Swisstopo')
ylabel('ALPNET')
xlim([-1.8 2.5])
ylim([-1.8 2.5])

%
% print(fig1, 'ALP_NETvsSwiss.eps','-depsc','-r300');
% print(fig1, 'ALP_NETvsSwiss.pdf','-dpdf','-r300');










