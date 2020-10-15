% cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP/Project_TippingPoints/VAR')
% save -v7.3 Fig1GC_Feldman.mat  diff_dT diff_sm diff_VOD diff_vpd diff_DSSF ...
%      dX_Dry_Mat SMAPCenterLongitudes SMAPCenterLatitudes cmapRedYelBu

%%
clc
clear

% cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP/Project_TippingPoints/VAR')
cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP/Project_TippingPoints/VAR/Figures/FinalFigures_GithubRepository')
load('Fig1GC_Feldman')
load('Fig1dXTimeEvolutionV2')

load coast   
clat = lat;                                                        
clon = long;

figure

dX_Dry_Mat1 = (dX_Dry_Mat./3).*7;
diff_DSSF1 = (diff_DSSF./3).*7;
diff_dT1 = (diff_dT./3).*7;
diff_sm1 = (diff_sm./3).*7;
diff_VOD1 = (diff_VOD./3).*7;
diff_vpd1 = (diff_vpd./3).*7;

% dX_Dry_Mat1 = dX_Dry_Mat;

for i = [3 6 9 12 15]
% for i = [6]

if i == 6
iplotA3 = subplot(5,3,i);  
% pcolor(SMAPCenterLongitudes,SMAPCenterLatitudes,...
%     dX_Dry_Mat(:,:,1)./diff_VOD(:,:,3)) ; shading flat
pcolor(SMAPCenterLongitudes,SMAPCenterLatitudes,...
    dX_Dry_Mat1(:,:,1)) ; shading flat
text(-29,-25,'$\frac{{\Delta}VOD}{{\Delta}t}$','Interpreter','latex','fontsize',22)
caxis([-0.07 0.07])
% caxis([-0.02 0.02])
elseif i == 3
iplotB3 = subplot(5,3,i);      
% pcolor(SMAPCenterLongitudes,SMAPCenterLatitudes,...
%     dX_Dry_Mat(:,:,2)./diff_sm(:,:,3)) ; shading flat
pcolor(SMAPCenterLongitudes,SMAPCenterLatitudes,...
    dX_Dry_Mat1(:,:,2)) ; shading flat
text(-22,-26,'$\frac{{\Delta}{\theta}}{{\Delta}t}$','Interpreter','latex','fontsize',22)
text(-25,-7,[{'Water-'};{'Limited'}],'fontsize',14)
caxis([-0.1 0.1])
% caxis([-0.02 0.02])
elseif i == 9
iplotC3 = subplot(5,3,i);  
% pcolor(SMAPCenterLongitudes,SMAPCenterLatitudes,...
% dX_Dry_Mat(:,:,3)./diff_DSSF(:,:,3)) ; shading flat
pcolor(SMAPCenterLongitudes,SMAPCenterLatitudes,...
dX_Dry_Mat1(:,:,3)) ; shading flat
text(-25,-25,'$\frac{{\Delta}Rs}{{\Delta}t}$','Interpreter','latex','fontsize',22)
caxis([-20 20])
% caxis([-5 5])
elseif i == 12
iplotD3 = subplot(5,3,i);      
% pcolor(SMAPCenterLongitudes,SMAPCenterLatitudes,...
%     dX_Dry_Mat(:,:,4)./diff_dT(:,:,3)) ; shading flat
pcolor(SMAPCenterLongitudes,SMAPCenterLatitudes,...
    dX_Dry_Mat1(:,:,4)) ; shading flat
text(-25,-25,'$\frac{{\Delta}dT}{{\Delta}t}$','Interpreter','latex','fontsize',22)
% caxis([-2 2])
caxis([-5 5])
elseif i == 15
iplotE3 = subplot(5,3,i);      
% pcolor(SMAPCenterLongitudes,SMAPCenterLatitudes,...
% dX_Dry_Mat(:,:,5)./diff_vpd(:,:,3)) ; shading flat
pcolor(SMAPCenterLongitudes,SMAPCenterLatitudes,...
dX_Dry_Mat1(:,:,5)) ; shading flat
text(-29,-25,'$\frac{{\Delta}VPD}{{\Delta}t}$','Interpreter','latex','fontsize',22)
% caxis([-0.5 0.5])
caxis([-1 1])
end

hold on
geoshow(clat,clon,'LineWidth',1,'Color','k')
ylim([-40 40])
xlim([-30 60])
% colormap(flipud(cmapRedYelBu))
colormap(cmapRedYelBu)
set(gca,'fontsize',14)
colorbar
set(gcf,'color','white')
set(gca,'yticklabel','')
set(gca,'xticklabel','')
% set(gca,'color',[0.9 0.9 0.9])
set(gca,'color',[0.8 0.8 0.8])

% caxis([-0.2 0.2])
end

for i = [2 5 8 11 14]
    
if i == 5
iplotA2 = subplot(5,3,i);
% dsmW = diff_VOD(:,:,1)./diff_VOD(:,:,3);
% dsmE = diff_VOD(:,:,2)./diff_VOD(:,:,3);
dsmW = diff_VOD1(:,:,1);
dsmE = diff_VOD1(:,:,2);
elseif i == 2
iplotB2 = subplot(5,3,i);
% dsmW = diff_sm(:,:,1)./diff_sm(:,:,3);
% dsmE = diff_sm(:,:,2)./diff_sm(:,:,3);
dsmW = diff_sm1(:,:,1);
dsmE = diff_sm1(:,:,2);
elseif i == 8
iplotC2 = subplot(5,3,i);
% dsmW = diff_DSSF(:,:,1)./diff_DSSF(:,:,3);
% dsmE = diff_DSSF(:,:,2)./diff_DSSF(:,:,3);
dsmW = diff_DSSF1(:,:,1);
dsmE = diff_DSSF1(:,:,2);
elseif i == 11
iplotD2 = subplot(5,3,i);
% dsmW = diff_dT(:,:,1)./diff_dT(:,:,3);
% dsmE = diff_dT(:,:,2)./diff_dT(:,:,3);
dsmW = diff_dT1(:,:,1);
dsmE = diff_dT1(:,:,2);
elseif i == 14
iplotE2 = subplot(5,3,i);    
% dsmW = diff_vpd(:,:,1)./diff_vpd(:,:,3);
% dsmE = diff_vpd(:,:,2)./diff_vpd(:,:,3);  
dsmW = diff_vpd1(:,:,1);
dsmE = diff_vpd1(:,:,2);  
end

missing = isnan(dsmW)|isnan(dsmE);
dsmW(missing)=nan; dsmE(missing)=nan;

% boxplot([dsmW(:),dsmE(:)],'colors',[0.8 0.8 0.8;0.2 0.8 0.2])
% boxplot([dsmW(:),dsmE(:)],'colors',[0 0 0])
boxplot([dsmE(:),dsmW(:)],'colors',[0 0 0])

set(findobj(gca,'type','line'),'linew',2)
h = findobj(gca,'Tag','Box');
for j=1:length(h)
  if j == 1
  p1 = patch(get(h(j),'XData'),get(h(j),'YData'),[0.3 1 0.3],'FaceAlpha',0.2);
  elseif j == 2
  p2 = patch(get(h(j),'XData'),get(h(j),'YData'),[0.6 0.6 0.6],'FaceAlpha',0.2);     
  end
end
h=findobj(gca,'tag','Outliers');
delete(h) 
% ylim([-0.2 0.2])

if i == 5
% ylim([-0.8 0.2])
end

if i == 5
% title('$\frac{{\Delta}VOD}{{\Delta}t}$','Interpreter','latex')
% text(1.15,0,'$\frac{{\Delta}VOD}{{\Delta}t}$','Interpreter','latex','fontsize',25)
% ylabel('Normalized Change Rate')
ylabel([{'\DeltaVOD/\Deltat'}; {'(Week^{-1})'}])
% ylim([-0.02 0.02])
ylim([-0.05 0.05])
% ylabel('\DeltaVOD/\Deltat')
elseif i == 2
ylabel([{'\Delta\theta/\Deltat'}; {'(m^3 m^{-3} Week^{-1})'}])
% text(1.25,0,'$\frac{{\Delta}{\theta}}{{\Delta}t}$','Interpreter','latex','fontsize',25)
% ylim([-0.02 0.02])
ylim([-0.12 0.12])
elseif i == 8
ylabel([{'\DeltaRs/\Deltat'} ;{'(W m^{-2} Week^{-1})'}])  
% text(1.25,0,'$\frac{{\Delta}Rs}{{\Delta}t}$','Interpreter','latex','fontsize',25)
% ylim([-5 5])
ylim([-40 40])
elseif i == 11
ylabel([{'\DeltadT/\Deltat'} ;{'(K Week^{-1})'}])  
% text(1.25,0,'$\frac{{\Delta}dT}{{\Delta}t}$','Interpreter','latex','fontsize',25)
% ylim([-1 1])
ylim([-5 5])
elseif i == 14
ylabel([{'\DeltaVPD/\Deltat'} ;{'(kPa Week^{-1})'}])  
% text(1.15,0,'$\frac{{\Delta}VPD}{{\Delta}t}$','Interpreter','latex','fontsize',25)
% ylim([-0.1 0.1])
ylim([-1 1])
end    

grid on
set(gca,'fontsize',14)
set(gca,'xticklabel','')
hold on 
plot([-2 5],[0 0],'--k','linewidth',2)
if i == 2
legend([p1 p2],'Water-Limited','Energy-Limited','location','north')
end

end



tthetaStar = find(SMQuantTimeEvo(2,:)<thetaStar); tthetaStar = min(tthetaStar);

tpulse = 0:3:12;
tpulseW = 0:3:(tthetaStar-1)*3;
tpulseE = ((tthetaStar-1)*3):3:12;
tpulseWInd = 1:1:(tthetaStar);
tpulseEInd = tthetaStar:1:5;


for i = [1 4 7 10 13]
    
if i == 4
timeMAT = VODQuantTimeEvo;
iplotA1 = subplot(5,3,i);
elseif i == 1
timeMAT = SMQuantTimeEvo;
iplotB1 = subplot(5,3,i);
elseif i == 7
timeMAT = RsQuantTimeEvo;
iplotC1 = subplot(5,3,i);
elseif i == 10
timeMAT = dTQuantTimeEvo;
iplotD1 = subplot(5,3,i);
elseif i == 13
timeMAT = VPDQuantTimeEvo;
iplotE1 = subplot(5,3,i);
end

patch([tpulseW fliplr(tpulseW)],...
    [timeMAT(1,tpulseWInd) fliplr(timeMAT(3,tpulseWInd))],[0.6 0.6 0.6],'facealpha',0.3)
hold on
patch([tpulseE fliplr(tpulseE)],...
    [timeMAT(1,tpulseEInd) fliplr(timeMAT(3,tpulseEInd))],[0.3 1 0.3],'facealpha',0.3)
plot(tpulse,timeMAT(2,:)','-k','linewidth',2)
hold on
plot([tpulse(tthetaStar) tpulse(tthetaStar)],[0 400],'--k','linewidth',1)
grid on
% plot([tthetaStar-1 tthetaStar-1],[nanmin(timeMAT(:)) nanmax(timeMAT(:))],'--k','linewidth',1)
% ylim([nanmin(timeMAT(:)) nanmax(timeMAT(:))])

if i == 4
ylabel('VOD')
% title('Plant Water Content','fontsize',27)
text(-7.5,0.2,[{'Plant'}; {'Water Content'}],'Rotation',90,'fontsize',16,...
    'horizontalalignment','center','fontweight','bold')
ylim([0.1 0.3])
elseif i == 1
ylabel('{\theta} (m^3 m^{-3})')
% text(5,0.03,[{'Water-Limited'}],'fontsize',14)
% text(0.1,0.04,[{'Energy-'};{'Limited'}],'fontsize',14)

text(3.5,0.235,[{'Water-Limited'}],'fontsize',14)
text(-1,0.25,[{'Energy-'};{'Limited'}],'fontsize',14)

text(3.2,0.04,[{'\theta^*'}],'fontsize',16)
text(-7.5,0.11,[{'Soil'}; {'Moisture'}],'Rotation',90,'fontsize',16,...
    'horizontalalignment','center','fontweight','bold')
% title('Soil Moisture','fontsize',27)
ylim([0.02 0.2])
elseif i == 7
ylabel('Rs (W m^{-2})')
% title('Solar Radiation','fontsize',27)
text(-7.5,250,[{'Solar'}; {'Radiation'}],'Rotation',90,'fontsize',16,...
    'horizontalalignment','center','fontweight','bold')
ylim([150 350])
elseif i == 10
ylabel('dT (K)')
% title([{'Diurnal Temperature'}; {'Amplitude'}],'fontsize',27)
text(-7.5,27.5,[{'Diurnal'};{'Temperature'};{'Amplitude'}],'Rotation',90,'fontsize',16,...
    'horizontalalignment','center','fontweight','bold')
ylim([15 40])
elseif i == 13
ylabel('VPD (kPa)')
% title([{'Vapor Pressure'}; {'Deficit'}],'fontsize',27)
text(-7.5,1.5,[{'Vapor Pressure'}; {'Deficit'}],'Rotation',90,'fontsize',16,...
    'horizontalalignment','center','fontweight','bold')
ylim([0.5 2.5])
xlabel('Time After Rainfall (Days)')
end

xlim([0 12])
set(gca,'fontsize',14)

end


% Boxplot (middle)
AMSA2 = get(iplotA2,'Position'); % top right fig
set(iplotA2,'Position',[AMSA2(1)-0.03 AMSA2(2)+0.02 AMSA2(3)+0.09 AMSA2(4)+0.02]);
AMSA2 = get(iplotA2,'Position'); % top right fig

AMSB2 = get(iplotB2,'Position'); % top right fig
set(iplotB2,'Position',[AMSA2(1) AMSB2(2)+0.04 AMSA2(3) AMSA2(4)]);
AMSB2 = get(iplotB2,'Position'); % top right fig

AMSC2 = get(iplotC2,'Position'); % top right fig
set(iplotC2,'Position',[AMSA2(1) AMSC2(2)+0.01 AMSA2(3) AMSA2(4)]);
AMSC2 = get(iplotC2,'Position'); % top right fig

AMSD2 = get(iplotD2,'Position'); % top right fig
set(iplotD2,'Position',[AMSA2(1) AMSD2(2)-0.01 AMSA2(3) AMSA2(4)]);
AMSD2 = get(iplotD2,'Position'); % top right fig

AMSE2 = get(iplotE2,'Position'); % top right fig
set(iplotE2,'Position',[AMSA2(1) AMSE2(2)-0.02 AMSA2(3) AMSA2(4)]);
AMSE2 = get(iplotE2,'Position'); % top right fig


% Time series (left)
AMSA1 = get(iplotA1,'Position'); % top right fig
set(iplotA1,'Position',[AMSA1(1)+0.03 AMSA2(2) AMSA1(3)-0.03 AMSA2(4)]);
AMSA1 = get(iplotA1,'Position'); % top right fig

AMSB1 = get(iplotB1,'Position'); % top right fig
set(iplotB1,'Position',[AMSA1(1) AMSB2(2) AMSA1(3) AMSA1(4)-0.02]);
AMSB1 = get(iplotB1,'Position'); % top right fig

AMSC1 = get(iplotC1,'Position'); % top right fig
set(iplotC1,'Position',[AMSA1(1) AMSC2(2) AMSA1(3) AMSA1(4)-0.005]);
AMSC1 = get(iplotC1,'Position'); % top right fig

AMSD1 = get(iplotD1,'Position'); % top right fig
set(iplotD1,'Position',[AMSA1(1) AMSD2(2) AMSA1(3) AMSA1(4)-0.005]);
AMSD1 = get(iplotD1,'Position'); % top right fig

AMSE1 = get(iplotE1,'Position'); % top right fig
set(iplotE1,'Position',[AMSA1(1) AMSE2(2) AMSA1(3) AMSA1(4)]);
AMSE1 = get(iplotE1,'Position'); % top right fig


% Maps (right)
AMSA3 = get(iplotA3,'Position'); % top right fig
set(iplotA3,'Position',[AMSA3(1)+0.045 AMSA2(2)-0.02 AMSA3(3)+0.02 AMSA2(4)+0.02]);
AMSA3 = get(iplotA3,'Position'); % top right fig

AMSB3 = get(iplotB3,'Position'); % top right fig
set(iplotB3,'Position',[AMSA3(1) AMSB2(2)-0.015 AMSA3(3) AMSA3(4)]);
AMSB3 = get(iplotB3,'Position'); % top right fig

AMSC3 = get(iplotC3,'Position'); % top right fig
set(iplotC3,'Position',[AMSA3(1) AMSC2(2)-0.025 AMSA3(3) AMSA3(4)]);
AMSC3 = get(iplotC3,'Position'); % top right fig

AMSD3 = get(iplotD3,'Position'); % top right fig
set(iplotD3,'Position',[AMSA3(1) AMSD2(2)-0.025 AMSA3(3) AMSA3(4)]);
AMSD3 = get(iplotD3,'Position'); % top right fig

AMSE3 = get(iplotE3,'Position'); % top right fig
set(iplotE3,'Position',[AMSA3(1) AMSE2(2)-0.03 AMSA3(3) AMSA3(4)]);
AMSE3 = get(iplotE3,'Position'); % top right fig

set(gcf,'position',[10 1 720 800])

% set(gcf,'position',[10 100 1400 750])

% annotation('line',[0 1],[0.79 0.79],'linewidth',3)
annotation('line',[0 1],[0.79 0.79],'linewidth',1.5,'linestyle','--')
annotation('line',[0 1],[0.605 0.605],'linewidth',1.5,'linestyle','--')
annotation('line',[0 1],[0.415 0.415],'linewidth',1.5,'linestyle','--')
annotation('line',[0 1],[0.225 0.225],'linewidth',1.5,'linestyle','--')
box on

annotation('line',[0 1],[0 0],'linewidth',3)
annotation('line',[0 0],[0 1],'linewidth',3) 
annotation('line',[1 1],[0 1],'linewidth',3)
annotation('line',[0 1],[1 1],'linewidth',3)

annotation('textarrow',[0.28 0.345],[0.95 0.95])
annotation('textarrow',[0.28 0.205],[0.95 0.95])
annotation('textarrow',[0.17 0.20],[0.95 0.95])
annotation('textarrow',[0.17 0.16],[0.95 0.95])

