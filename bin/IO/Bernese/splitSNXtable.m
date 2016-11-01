function [Table] = splitSNXtable(SNXTable)
    SNXTable = SolutionEstimateTable;
    
    StationList = unique(SNXTable(:,3));
    NumberOfStations = length(StationList);
    
    Table = zeros((size(SNXTable,1)/6),12);
    
    for block = 1:size(SNXTable,1)/6
    row = block*6-6+1;
        Station = SNXTable(row,3);
        PT = SNXTable(row,4);
        Sol = cell2mat(SNXTable(row,5));
        CRD = cell2mat(SNXTable(row:row+2,6));
        SigmaCRD = cell2mat(SNXTable(row:row+2,7));
        VEL = cell2mat(SNXTable(row+3:row+5,6));
        SigmaVEL = cell2mat(SNXTable(row+3:row+5,7));
        
        Table(block,:) = [CRD',VEL', SigmaCRD', SigmaVEL'];
        
    end

    %%
    
    [V_vertical] = VerticalProjection(Table(:,1),Table(:,2),Table(:,3),Table(:,4),Table(:,5),Table(:,6));
    clr = lines(5);
    counterS = 0;
    counterU = 0;
    for i=1:size(V_vertical)
       V_mag(i) = norm(V_vertical(i,:));
       e_r   = Table(i,1:3)./norm(Table(i,1:3));
       e_Vup = V_vertical(i,:)./norm(V_vertical(i,:));
       if norm(e_r + e_Vup) > 1  % && V_mag(i) < 0.004    
           counterU = counterU + 1 ;
           Flag_U(counterU) = i;
       elseif norm(e_r + e_Vup) < 1  % && V_mag(i) < 0.004
           counterS = counterS + 1 ;
           Flag_S(counterS) = i; 
       end
    end
    
    Flag_all = sort([Flag_S, Flag_U]);

    %%
    close all
    figure(1)
    hold on
    axis equal
    grid on
    s = 1000*1000*50;
    plot3(Table(:,1),Table(:,2),Table(:,3),'.b')
%     quiver3(Table(:,1),Table(:,2),Table(:,3),Table(:,4)*s,Table(:,5)*s,Table(:,6)*s,0,'r')
    quiver3(Table(Flag_S,1),Table(Flag_S,2),Table(Flag_S,3),V_vertical(Flag_S,1)*s,V_vertical(Flag_S,2)*s,V_vertical(Flag_S,3)*s,0,'Color',clr(3,:),'LineWidth', 2)
    quiver3(Table(Flag_U,1),Table(Flag_U,2),Table(Flag_U,3),V_vertical(Flag_U,1)*s,V_vertical(Flag_U,2)*s,V_vertical(Flag_U,3)*s,0,'Color',clr(2,:),'LineWidth', 2)
    
    title('Velocity field')
    xlabel('X ,[m]')
    ylabel('Y ,[m]')
    zlabel('Z ,[m]')
    
%     close all
%     figure(2)
%     hold on
%     axis equal
%     grid on
    for i = 1:size(Table,1)
        SigmaCRD = Table(i,7:9)*10^1;
        SigmaVEL = Table(i,10:12)*10^9;
          [x ,y ,z] = ellipsoid(Table(i,1),Table(i,2),Table(i,3),SigmaCRD(1),SigmaCRD(2),SigmaCRD(3),12);
          surf(x ,y ,z, ones(size(x)))
          [x ,y ,z] = ellipsoid(Table(i,1),Table(i,2),Table(i,3),SigmaVEL(1),SigmaVEL(2),SigmaVEL(3),12);
          surf(x ,y ,z, 2*ones(size(x)))        
          colormap([0.5  0.5  0.5; 1 0 0])

    
    end
    %%
    close all
    clc
   [x ,y ,z] = ellipsoid(Table(i,1),Table(i,2),Table(i,3),SigmaCRD(1),SigmaCRD(2),SigmaCRD(3),12);
   surf(x ,y ,z, 0.5*zeros(size(x)))
   colormap([0.5  0.5  0.5])


    %%

end