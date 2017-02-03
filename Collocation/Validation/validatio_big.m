%% load and prepare swisstopo
validation_swiss

%% load EPN
load('../../dat/SNX/EPN_A_IGb08_no_COVA.SNX.mat')

%% France

load('Validation/FranceVELtable.mat')
FRvel = FranceVELtable;

names_comm_AlFr = intersect(names, FRvel.Sites);




%% remove plare motion from SNX solution

% Average vel

[V_res_snx, V_pl_snx] = remove_plate_motion(EPN_SNX.SOLUTION.ESTIMATE.Data.CRD, EPN_SNX.SOLUTION.ESTIMATE.Data.VEL, Omega_Eur);
[Ve_res_epn, Vn_res_epn, Vu_res_epn, lat_epn,lon_epn] = XYZ2ENU(EPN_SNX.SOLUTION.ESTIMATE.Data.CRD,V_res_snx);
[Ve_pl_epn,  Vn_pl_epn,  Vu_pl_epn ]              = XYZ2ENU(EPN_SNX.SOLUTION.ESTIMATE.Data.CRD,V_pl_snx);

names_epn = EPN_SNX.SOLUTION.ESTIMATE.StationData.Site;
SolN_epn  = EPN_SNX.SOLUTION.ESTIMATE.StationData.SolN;
range_s = 1:length(names_epn);
iiEPNcom = range_s(ismember(names_epn, names));
iiALPcom = range_s(ismember(names, names_epn));

% table_EPN = table(lat_s,lon_s, Ve_res_s, Vn_res_s, Vu_res_s, names_epn, 'RowNames', [names_epn, SolN_epn])

%%  compute difference


v_epn_mean =  zeros(length(iiALPcom),5);
range_c = 1:length(names_epn);
for i = 1:length(iiALPcom)
   iiCom = range_c(ismember( names_epn, names(iiALPcom(i)) ) );
    names_epn_uniq(i,1) = names(iiALPcom(i));
    v_epn_mean(i,:) = [long(iiALPcom(i)), lat(iiALPcom(i)), mean(Ve_res_epn(iiCom)), mean(Vn_res_epn(iiCom)), mean(Vu_res_epn(iiCom))]; 
end

names_comm_AlEp = names(iiALPcom);
table_Alp2 = table(long, lat, Ve_res*1000, Vn_res*1000, Vu_res*1000, names, 'RowNames', names)
table_Alp2.Properties.VariableNames = {'Long', 'Lat', 'Ve_res', 'Vn_res', 'Vu', 'Site'};
table_EPN1 = table(v_epn_mean(:,1),v_epn_mean(:,2),v_epn_mean(:,3)*1000,v_epn_mean(:,4)*1000,v_epn_mean(:,5)*1000,names_epn_uniq, 'RowNames', names_epn_uniq)
table_EPN1.Properties.VariableNames = {'Long', 'Lat', 'Ve_res', 'Vn_res', 'Vu', 'Site'};

%% trep alp table Vres

alp_table = table(lat, long, Ve*1000, Vn*1000, Vu*1000, Ve_res*1000, Vn_res*1000, Vu_res*1000, names_alp, 'RowNames', names_alp)



%%
    names_epn_51 = table_EPN1{:,6}

    names_comm_es = intersect(names_common,names_epn_51)
    
%%

% Defaults for this blog post
width =  7;     % Width in inches
height = 7;    % Height in inches
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
fig2 = figure(2);

% ALPNET vs swisstopo
subplot(3,3,1)
hold on;  axis equal
plot([-1 1.5],[-1 1.5],'-b')
plot( alp_table{names_common,6} , table_sw_res{:,3}, 'or','MarkerFaceColor','r') 
text(-.7, 1.2,   [ 'rms = ',  num2str(rms( alp_table{names_common,6} - table_sw_res{:,3}), '%5.2f'), ' mm/yr'])
title(['V_E_a_s_t, mm/yr'])
ylabel('Swisstopo')
xlabel('ALPNET')
xlim([-1 1.5])
ylim([-1 1.5])

subplot(3,3,2)
hold on; axis equal
plot([-1 4],[-1 4],'-b')
plot( alp_table{names_common,7} , table_sw_res{:,4} , 'or','MarkerFaceColor','r') 
text(-0.8, 3.2,   [ 'rms = ',  num2str(rms( alp_table{names_common,7} - table_sw_res{:,4} ), '%5.2f'), ' mm/yr'])
title(['V_N_o_r_t_h, mm/yr'])
ylabel('Swisstopo')
xlabel('ALPNET')
xlim([-1 3.5])
ylim([-1 3.5])

subplot(3,3,3)
hold on;  axis equal
plot([-1.8 2.5],[-1.8 2.5],'-b')
plot( alp_table{names_common,5} , table_sw_res{:,5} , 'or', 'MarkerFaceColor','r') 
text(-1.6, 2.1, [ 'rms = ',  num2str(rms( alp_table{names_common,5} - table_sw_res{:,5} ), '%5.2f'), ' mm/yr'])
title(['V_U_p, mm/yr'])
ylabel('Swisstopo')
xlabel('ALPNET')
xlim([-1.8 2.5])
ylim([-1.8 2.5])


% ALPNET vs EPN
subplot(3,3,4)
hold on;  axis equal
plot([-1 1.5],[-1 1.5],'-b')
plot( table_Alp2{names_comm_AlEp,3} , table_EPN1{:,3} , 'or','MarkerFaceColor','r') 
text(-0.7, 1.3,   [ 'rms = ',  num2str(rms( table_Alp2{names_comm_AlEp,3} - table_EPN1{:,3} ), '%5.2f'), ' mm/yr'])
title(['V_E_a_s_t, mm/yr'])
ylabel('EPN')
xlabel('ALPNET')
xlim([-1 1.5])
ylim([-1 1.5])

subplot(3,3,5)
hold on; axis equal
plot([-1 4],[-1 4],'-b')
plot( table_Alp2{names_comm_AlEp,4} , table_EPN1{:,4} , 'or','MarkerFaceColor','r') 
text(-0.7, 3.2,   [ 'rms = ',  num2str(rms( table_Alp2{names_comm_AlEp,4} - table_EPN1{:,4}), '%5.2f'), ' mm/yr'])
title(['V_N_o_r_t_h, mm/yr'])
ylabel('EPN')
xlabel('ALPNET')
xlim([-1 3.5])
ylim([-1 3.5])

subplot(3,3,6)
hold on;  axis equal
plot([-1.8 2.5],[-1.8 2.5],'-b')
plot(  table_Alp2{names_comm_AlEp,5}, table_EPN1{:,5}, 'or', 'MarkerFaceColor','r') 
text(-1.6, 2.1, [ 'rms = ',  num2str(rms( table_EPN1{:,5} - table_Alp2{names_comm_AlEp,5} ), '%5.2f'), ' mm/yr'])
title(['V_U_p, mm/yr'])
ylabel('EPN')
xlabel('ALPNET')
xlim([-1.8 2.5])
ylim([-1.8 2.5])

% ALPNET vs France
subplot(3,3,7)
hold on;  axis equal
plot([-1 1.5],[-1 1.5],'-b')
plot( table_Alp2{names_comm_AlFr,3} ,FRvel{names_comm_AlFr,4} , 'or','MarkerFaceColor','r') 
text(-0.7, 1.3,   [ 'rms = ',  num2str(rms( FRvel{names_comm_AlFr,4} - table_Alp2{names_comm_AlFr,3}), '%5.2f'), ' mm/yr'])
title(['V_E_a_s_t, mm/yr'])
ylabel('France')
xlabel('ALPNET')
xlim([-1 1.5])
ylim([-1 1.5])


subplot(3,3,8)
hold on; axis equal
plot([-1 4],[-1 4],'-b')
plot( table_Alp2{names_comm_AlFr,4}, FRvel{names_comm_AlFr,5}, 'or','MarkerFaceColor','r') 
text(-0.7, 3.2,   [ 'rms = ',  num2str(rms( FRvel{names_comm_AlFr,5} -  table_Alp2{names_comm_AlFr,4}), '%5.2f'), ' mm/yr'])
title(['V_N_o_r_t_h, mm/yr'])
ylabel('France')
xlabel('ALPNET')
xlim([-1 3.5])
ylim([-1 3.5])

subplot(3,3,9)
hold on;  axis equal
plot([-1.8 2.5],[-1.8 2.5],'-b')
plot( table_Alp2{names_comm_AlFr,5}, FRvel{names_comm_AlFr,6} , 'or','MarkerFaceColor','r') 
text(-1.6, 2.1, [ 'rms = ',  num2str(rms( FRvel{names_comm_AlFr,6} - table_Alp2{names_comm_AlFr,5}), '%5.2f'), ' mm/yr'])
title(['V_U_p, mm/yr'])
ylabel('France')
xlabel('ALPNET')
xlim([-1.8 2.5])
ylim([-1.8 2.5])

%%
print(fig2, '-depsc','-r300', 'CompareVelAll.eps')
print(fig2, '-dpdf','-r300', 'CompareVelAll.pdf')





%% EPN wv swisstopo

% % ALPNET vs EPN
% subplot(3,3,7)
% hold on;  axis equal
% plot([-1 1.5],[-1 1.5],'-b')
% Ve_epn_comm = table_EPN1{names_comm_es,3};
% Ve_alp2_comm = table_sw_res{names_comm_es,3};
% plot( Ve_epn_comm , Ve_alp2_comm , 'or','MarkerFaceColor','r') 
% title(['V_E_a_s_t, mm/yr'])
% % text(17, 21.5, [ 'mean(dV_e) = ', num2str(mean(Ve_ch_comm - Ve_al_comm), '%5.3f'), ' mm/yr'])
% text(-0.7, 1.3,   [ 'rms = ',  num2str(rms( Ve_epn_comm - Ve_alp2_comm), '%5.2f'), ' mm/yr'])
% xlabel('EPN')
% ylabel('ALPNET')
% xlim([-1 1.5])
% ylim([-1 1.5])
% 
% 
% subplot(3,3,8)
% hold on; axis equal
% plot([-1 4],[-1 4],'-b')
% Vn_epn_comm = table_EPN1{names_comm_es,4};
% Vn_alp2_comm = table_sw_res{names_comm_es,4};
% plot( Vn_epn_comm , Vn_alp2_comm , 'or','MarkerFaceColor','r') 
% title(['V_N_o_r_t_h, mm/yr'])
% % text(15.2, 17.8, [ 'mean(dV_e) = ', num2str(mean(Vn_ch_comm - Vn_al_comm), '%5.3f'), ' mm/yr'])
% text(-0.7, 3.5,   [ 'rms = ',  num2str(rms( Vn_epn_comm - Vn_alp2_comm), '%5.2f'), ' mm/yr'])
% xlabel('EPN')
% ylabel('ALPNET')
% xlim([-1 4])
% ylim([-1 4])
% 
% subplot(3,3,9)
% hold on;  axis equal
% plot([-1.8 2.5],[-1.8 2.5],'-b')
% Vu_epn_comm = table_sw_res{names_comm_es,5};
% Vu_alp2_comm = table_Alp2{names_comm_es,5};
% plot( Vu_epn_comm , Vu_alp2_comm , 'or', 'MarkerFaceColor','r') 
% title(['V_U_p, mm/yr'])
% text(-1.6, 2.1, [ 'rms = ',  num2str(rms( Vu_epn_comm - Vu_alp2_comm), '%5.2f'), ' mm/yr'])
% xlabel('EPN')
% ylabel('swisstopo')
% xlim([-1.8 2.5])
% ylim([-1.8 2.5])




