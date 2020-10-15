
% cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP/Project_TippingPoints/VAR')
% save -v7.3 Fig2Part1GC_Feldman.mat beta_LS_Dry_Mat1
% save -v7.3 Fig2Part2GC_Feldman.mat beta_LS_Dry_Mat2
% save -v7.3 Fig2Part3GC_Feldman.mat p_GC_CV_Mat
% save -v7.3 Fig2Part4GC_Feldman.mat dX_Dry_Mat AfricaFlag ...
%      SMAPCenterLongitudes SMAPCenterLatitudes cmapRedYelBu MeanMAP
%  
% cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP Data/General Data')
% load('GPMMean_201504_201903')
% beta_LS_Dry_Mat1 = beta_LS_Dry_Mat(:,:,1:12);
% beta_LS_Dry_Mat2 = beta_LS_Dry_Mat(:,:,13:25);

%% VOD Betas (GC)

clc
clear
%%%% Change this repository to where you have the files stored! %%%%
cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP/Project_TippingPoints/VAR/Figures/FinalFigures_GithubRepository')
% cd('/Users/andrewfeldman/Dropbox(MIT)/SMAP/Project_TippingPoints/VAR')
load('Fig2Part1GC_Feldman')
load('Fig2Part2GC_Feldman')
load('Fig2Part3GC_Feldman')
load('Fig2Part4GC_Feldman')

beta_LS_Dry_Mat = cat(3,beta_LS_Dry_Mat1,beta_LS_Dry_Mat2);

load coast   
clat = lat;                                                        
clon = long; 

Row_Af = 300:1300;
Col_Af = 1700:2500;
MeanMAP = MeanMAP(Row_Af,Col_Af);

p_Het_Dry_Flag = p_GC_CV_Mat.*AfricaFlag;
p_Het_Dry_Flag(p_Het_Dry_Flag>0.05)=nan;
p_Het_Dry_Flag(isfinite(p_Het_Dry_Flag))=1;

p_Het_Dry_NonSig = p_GC_CV_Mat.*AfricaFlag;
p_Het_Dry_NonSig(p_Het_Dry_NonSig<0.05)=nan;
p_Het_Dry_NonSig(isfinite(p_Het_Dry_NonSig))=1;

HistMat_Dry = beta_LS_Dry_Mat.*p_Het_Dry_Flag;
HistFin_Dry = HistMat_Dry; HistFin_Dry(isfinite(HistFin_Dry))=1;

nSamp_Dry = beta_LS_Dry_Mat(:,:,1).*AfricaFlag;
nSamp_Dry(isfinite(nSamp_Dry))=1;
PerBetaSig_Dry = nan(25,1);
for i = 1:length(PerBetaSig_Dry)
    A = HistFin_Dry(:,:,i);
    PerBetaSig_Dry(i) = nansum(A(:))./nansum(nSamp_Dry(:));
end

dX_Dry_Mat1 = dX_Dry_Mat.*(7/3);

figure
%%%%% Drydowns Only WITH VPD %%%%%%
for iind = 2:5
betaMatdXSig = beta_LS_Dry_Mat(:,:,iind).*dX_Dry_Mat1(:,:,iind);
if iind == 2
iplotA = subplot(2,2,iind-1);
elseif iind == 3
iplotB = subplot(2,2,iind-1);
elseif iind == 4
iplotC = subplot(2,2,iind-1);
elseif iind == 5
iplotD = subplot(2,2,iind-1);
end

pcolor(SMAPCenterLongitudes,SMAPCenterLatitudes,...
    betaMatdXSig) ; shading flat
hold on
% cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP/Project_PlantSoilVODStudyExtended')
stipple(SMAPCenterLongitudes,SMAPCenterLatitudes,p_Het_Dry_Flag(:,:,iind)==1,'density',150,'markersize',4) ; shading flat

hold on
set(gca,'yticklabel','')
set(gca,'xticklabel','')

if iind == 2
% title('\theta_{t-1} \Rightarrow VOD_{t}')
title('${{\Delta}VOD_{{\theta_{t-1}}{\rightarrow}VOD_{t}}}$',...
    'Interpreter','latex')
colormap(iplotA,cmapRedYelBu)
text(39,-29,'\downarrow\theta_{t-1}','fontsize',16,'FontWeight','bold')
text(37,-36,'\downarrowVOD_{t}','fontsize',16,'FontWeight','bold')
elseif iind == 3
% title('Rs_{t-1} \Rightarrow VOD_{t}')
title('${{\Delta}VOD_{{Rs_{t-1}}{\rightarrow}VOD_{t}}}$',...
    'Interpreter','latex')
colormap(iplotB,cmapRedYelBu)
elseif iind == 4
% title('dT_{t-1} \Rightarrow VOD_{t}')
title('${{\Delta}VOD_{{dT_{t-1}}{\rightarrow}VOD_{t}}}$',...
    'Interpreter','latex')
colormap(iplotC,cmapRedYelBu)
text(38,-29,'\uparrowdT_{t-1}','fontsize',16,'FontWeight','bold')
text(37,-36,'\downarrowVOD_{t}','fontsize',16,'FontWeight','bold')
elseif iind == 5
% title('VPD_{t-1} \Rightarrow VOD_{t}')
title('${{\Delta}VOD_{{VPD_{t-1}}{\rightarrow}VOD_{t}}}$',...
    'Interpreter','latex')
colormap(iplotD,cmapRedYelBu)
end

geoshow(clat,clon,'LineWidth',1,'Color','k')
ylim([-40 40])
xlim([-30 60])
set(gca,'fontsize',18)
colorbar
set(gcf,'color','white')
% caxis([-0.01 0.01])
caxis([-0.025 0.025])
set(gca,'color',[0.75 0.75 0.75])
% set(gca,'color',[0.8 0.8 0.8])

if iind == 2
hiDiag = axes('Position',[.095 .575 .1 .13]);
elseif iind == 3
hiDiag = axes('Position',[.575 .575 .1 .13]);
elseif iind == 4
hiDiag = axes('Position',[.095 .09 .1 .13]);
elseif iind == 5
hiDiag = axes('Position',[.575 .09 .1 .13]);
end

ASig = beta_LS_Dry_Mat(:,:,iind).*p_Het_Dry_Flag(:,:,iind)...
    .*dX_Dry_Mat1(:,:,iind);
ANoSig = beta_LS_Dry_Mat(:,:,iind).*p_Het_Dry_NonSig(:,:,iind)...
    .*dX_Dry_Mat1(:,:,iind);

[k1, f1] = ksdensity(ASig(:),'bandwidth',0.00005);
% [k2, f2] = ksdensity(ANoSig(:),'bandwidth',0.0005);
plot(f1,k1,'linewidth',1.5,'color','k')
hold on
% plot(f2,k2,'linewidth',1.5,'color','k')
% hold on
patch([f1 fliplr(f1)],[k1 zeros(size(k1))],[0.5 0.5 0.5],'facealpha',0.4)
% hold on
% patch([f2 fliplr(f2)],[k2 zeros(size(k2))],[1 1 1],'facealpha',0.1)
set(gca,'yticklabel','')
grid on

% hist(ASig(:),100);
% hold on
% hist(ANoSig(:),100,'k');

% title(['Sig = ' sprintf('%0.0f',PerBetaSig_Dry(iind)*100) '%'])
% title([{[sprintf('%0.0f',PerBetaSig_Dry(iind)*100) '% Significant']}])
% title([{[sprintf('%0.0f',PerBetaSig_Dry(iind)*100) '% (p<0.05)']}])
title([{[sprintf('%0.0f',PerBetaSig_Dry(iind)*100) '%']}])
% xlim([-0.01 0.01])
xlim([-0.025 0.025])
colormap(hiDiag,gray)
set(gca,'fontsize',12)

end

AMSA = get(iplotA,'Position'); % top right fig
set(iplotA,'Position',[AMSA(1)-0.06 AMSA(2)-0.04 AMSA(3)+0.06 AMSA(4)+0.06]);
AMSA = get(iplotA,'Position'); % top right fig

AMSB = get(iplotB,'Position'); % top right fig
set(iplotB,'Position',[AMSB(1)-0.02 AMSA(2) AMSA(3) AMSA(4)]);
AMSB = get(iplotB,'Position'); % top right fig

AMSC = get(iplotC,'Position'); % top right fig
set(iplotC,'Position',[AMSA(1) AMSC(2)-0.05 AMSA(3) AMSA(4)]);
AMSC = get(iplotC,'Position'); % top right fig

AMSD = get(iplotD,'Position'); % top right fig
set(iplotD,'Position',[AMSB(1) AMSC(2) AMSA(3) AMSA(4)]);
AMSD = get(iplotD,'Position'); % top right fig

Atext = ['A'];
annotation('textbox',[0.02 0.9 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',24,'EdgeColor','none')
Atext = ['B'];
annotation('textbox',[0.49 0.9 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',24,'EdgeColor','none')
Atext = ['C'];
annotation('textbox',[0.02 0.42 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',24,'EdgeColor','none')
Atext = ['D'];
annotation('textbox',[0.49 0.42 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',24,'EdgeColor','none')

set(gcf,'position',[100 100 800 650])

annotation('line',[0 1],[0 0],'linewidth',3)
annotation('line',[0 0],[0 1],'linewidth',3)
annotation('line',[1 1],[0 1],'linewidth',3)
annotation('line',[0 1],[1 1],'linewidth',3)


%%%%%% Interactions Figure 3 %%%%%%

AfricaAll = beta_LS_Dry_Mat(:,:,1).*AfricaFlag;
AfricaAll(isfinite(AfricaAll))=1;

% pFlag = p_NW_Dry_Mat;
% pFlag = p_Hom_Dry_Mat;
% pFlag = p_Het_Dry_Mat;
pFlag = p_GC_CV_Mat;
% pFlag = p_Boot_Dry_Mat;
% pFlag = p_BlockBoot_Dry_Mat;
pFlag(pFlag>0.05)=nan; pFlag(isfinite(pFlag))=1;

% beta_LS_Dry_Mat_Sig = beta_LS_Blockboot_Dry.*pFlag;
% beta_LS_Dry_Mat_Sig = beta_LS_boot_Dry.*pFlag;
beta_LS_Dry_Mat_Sig = beta_LS_Dry_Mat.*pFlag;

%dVOD
VOD_X = beta_LS_Dry_Mat_Sig(:,:,1:5).*dX_Dry_Mat.*AfricaFlag; 
SM_X = beta_LS_Dry_Mat_Sig(:,:,6:10).*dX_Dry_Mat.*AfricaFlag;
Rs_X = beta_LS_Dry_Mat_Sig(:,:,11:15).*dX_Dry_Mat.*AfricaFlag;
dT_X = beta_LS_Dry_Mat_Sig(:,:,16:20).*dX_Dry_Mat.*AfricaFlag;
VPD_X = beta_LS_Dry_Mat_Sig(:,:,21:25).*dX_Dry_Mat.*AfricaFlag;

VOD_X1 =VOD_X; % save vars
% VOD_X1(VOD_X1>0)=nan;
% VOD_X1(VOD_X1<0)=nan;

SMPos = VOD_X1(:,:,2);
SMPos(SMPos<0)=nan; SMPos(isfinite(SMPos))=1;
SMNeg = VOD_X1(:,:,2);
SMNeg(SMNeg>0)=nan; SMNeg(isfinite(SMNeg))=1;
PerNegSM = nansum(SMNeg(:))./(nansum(SMPos(:))+nansum(SMNeg(:)));

dTPos = VOD_X1(:,:,4);
dTPos(dTPos<0)=nan; dTPos(isfinite(dTPos))=1;
dTNeg = VOD_X1(:,:,2);
dTNeg(dTNeg>0)=nan; dTNeg(isfinite(dTNeg))=1;
PerNegdT = nansum(dTNeg(:))./(nansum(dTPos(:))+nansum(dTNeg(:)));

% Number of sig + VOD<0
VOD_X(isfinite(VOD_X))=1; % All VOD %
% VOD_X(VOD_X>0)=nan; VOD_X(isfinite(VOD_X))=1; % Decrease VOD %
% VOD_X(VOD_X<0)=nan; VOD_X(isfinite(VOD_X))=1; % Increase VOD %
SM_X(SM_X>0)=nan; SM_X(isfinite(SM_X))=1; % Decrease SM
Rs_X(Rs_X<0)=nan; Rs_X(isfinite(Rs_X))=1; % Increase Rs
dT_X(dT_X<0)=nan; dT_X(isfinite(dT_X))=1; % Increase dT
VPD_X(VPD_X<0)=nan; VPD_X(isfinite(VPD_X))=1; % Increase VPD

%Compute Interactions
VOD_X_Int = beta_LS_Dry_Mat_Sig(:,:,1:5).*AfricaFlag;
SM_X_Int = beta_LS_Dry_Mat_Sig(:,:,6:10).*dX_Dry_Mat.*AfricaFlag;
dT_X_Int = beta_LS_Dry_Mat_Sig(:,:,16:20).*dX_Dry_Mat.*AfricaFlag;
Rs_X_Int = beta_LS_Dry_Mat_Sig(:,:,11:15).*dX_Dry_Mat.*AfricaFlag;
VPD_X_Int = beta_LS_Dry_Mat_Sig(:,:,21:25).*dX_Dry_Mat.*AfricaFlag;

VOD_SM_XInt = VOD_X_Int(:,:,2).*SM_X_Int; %1 3 4 5
VOD_dT_XInt = VOD_X_Int(:,:,4).*dT_X_Int; %1 2 3 5
VOD_Rs_XInt = VOD_X_Int(:,:,3).*Rs_X_Int; %1 2 4 5
VOD_VPD_XInt = VOD_X_Int(:,:,5).*VPD_X_Int; %1 2 3 4
TotalInteractions = nansum(VOD_SM_XInt(:,:,[1 3 4 5]),3)+...
    nansum(VOD_dT_XInt(:,:,[1 2 3 5]),3)+...
    nansum(VOD_Rs_XInt(:,:,[1 2 4 5]),3)+...
    nansum(VOD_VPD_XInt(:,:,[1 2 3 4]),3);

% VOD_SM_XInt(VOD_SM_XInt>0)=nan; 
% VOD_SM_XInt(VOD_SM_XInt<0)=nan; 
VOD_SM_XInt1 = VOD_SM_XInt;
VOD_SM_XInt(isfinite(VOD_SM_XInt))=1; % Decrease VOD

% VOD_dT_XInt(VOD_dT_XInt>0)=nan; 
% VOD_dT_XInt(VOD_dT_XInt<0)=nan; 
VOD_dT_XInt1 = VOD_dT_XInt; % Save vars
VOD_dT_XInt(isfinite(VOD_dT_XInt))=1; % Decrease VOD

% VOD_SM_XInt(VOD_SM_XInt<0)=nan; VOD_SM_XInt(isfinite(VOD_SM_XInt))=1; % Decrease VOD
% VOD_dT_XInt(VOD_dT_XInt<0)=nan; VOD_dT_XInt(isfinite(VOD_dT_XInt))=1; % Decrease VOD

VOD_SM_Frac = VOD_SM_XInt1./VOD_X1(:,:,2); %Indirect/Direct SM->VOD
VOD_dT_Frac = VOD_dT_XInt1./VOD_X1(:,:,4); %Indirect/Direct dT->VOD

MapBound = 0:200:2000;
MapBoundAve = nanmean(cat(1,MapBound(1:end-1),MapBound(2:end)),1);

% VOD_SM_Frac = (VOD_SM_XInt1./VOD_X1(:,:,2)); %Indirect/Direct SM->VOD
% VOD_dT_Frac = (VOD_dT_XInt1./VOD_X1(:,:,4)); %Indirect/Direct dT->VOD

% MapBound = 0:250:2000;
MapBound = 0:200:2000;
pProbVOD_SigNeg = nan(4,length(MapBound)-1);
pProbSM_SigNeg = nan(4,length(MapBound)-1);
pProbdT_SigNeg = nan(4,length(MapBound)-1);

pProbVOD_SM_SigNeg = nan(5,length(MapBound)-1);
pProbVOD_dT_SigNeg = nan(5,length(MapBound)-1);

for iP = 1:5
    for iMAP = 1:length(MapBound)-1
    % Only selected region Sig
    ASigNegSM = VOD_SM_XInt(:,:,iP);
    ASigNegSM(MeanMAP<MapBound(iMAP)|MeanMAP>=MapBound(iMAP+1))=nan;

    ASigNegdT = VOD_dT_XInt(:,:,iP);
    ASigNegdT(MeanMAP<MapBound(iMAP)|MeanMAP>=MapBound(iMAP+1))=nan;

    % All available locations
    Atot = AfricaAll;
    Atot(MeanMAP<MapBound(iMAP)|MeanMAP>=MapBound(iMAP+1))=nan;

    pProbVOD_SM_SigNeg(iP,iMAP) = nansum(ASigNegSM(:))./nansum(Atot(:))*100;        
    pProbVOD_dT_SigNeg(iP,iMAP) = nansum(ASigNegdT(:))./nansum(Atot(:))*100;            
    end
end
  
% VOD
for iP = 1:4
    for iMAP = 1:length(MapBound)-1
    % Only selected region Sig
    ASigNeg = VOD_X(:,:,iP+1);
    ASigNeg(MeanMAP<MapBound(iMAP)|MeanMAP>=MapBound(iMAP+1))=nan;

    % All available locations
    Atot = AfricaAll;
    Atot(MeanMAP<MapBound(iMAP)|MeanMAP>=MapBound(iMAP+1))=nan;
    
    pProbVOD_SigNeg(iP,iMAP) = nansum(ASigNeg(:))./nansum(Atot(:))*100;    
    end
end

MapBoundAve = nanmean(cat(1,MapBound(1:end-1),MapBound(2:end)),1);

figure
iplotA = subplot(2,2,1);
plot(MapBoundAve,pProbVOD_SigNeg(1,:),'-o','color',[0 0 0],'linewidth',3)
hold on
colorvec = [0 0.6 0; 0 0 0; 1 0.7 0.1; 0.3 0.7 0.9; 1 0 0.1];
for i = [1 3:5]
plot(MapBoundAve,pProbVOD_SM_SigNeg(i,:),'-o','linewidth',2,'color',...
    colorvec(i,:))
hold on
end

legend('\theta_{t-1}\rightarrow VOD_{t} (Direct)',...
    'VOD_{t-2}\rightarrow \theta_{t-1}\rightarrow VOD_{t}',...
    '   Rs_{t-2}\rightarrow \theta_{t-1}\rightarrow VOD_{t}',...
    '   dT_{t-2}\rightarrow \theta_{t-1}\rightarrow VOD_{t}',...
    'VPD_{t-2}\rightarrow \theta_{t-1}\rightarrow VOD_{t}','fontsize',12,'position',[0.31 0.8 0.1 0.1])
% ylabel([{'Plant Drying Influence'}; {'via Soil Moisture (% Area)'}])
% ylabel([{'Significant VOD'}; {'Influence (% Area)'}])
% ylabel([{'Significant Influence'}; {'on VOD (% Area)'}])
% ylabel([{'Interaction'};{'Areal Prevalence (%)'}])
% ylabel([{'Frequency of Statistically'};{'Significant Interaction (% Area)'}])
ylabel([{'Prevalence of Statistically'};{'Significant Interaction (% Area)'}])

ylim([0 70])
grid on
set(gca,'fontsize',14)
xlabel('Mean Annual Precipitation (mm)')
% text(2200,90,[{'Influence of Soil Moisture Interactions on Plant Water Content'}],...
%     'fontsize',20,'horizontalalignment', 'center','FontWeight','bold')
text(2300,77,[{'Effect of Soil Moisture Interactions on Plant Water Content'}],...
    'fontsize',17,'horizontalalignment', 'center','FontWeight','bold')

% title('Soil Moisture \rightarrow VOD')

iplotC = subplot(2,2,3);
plot(MapBoundAve,pProbVOD_SigNeg(3,:),'-o','color',[0 0 0],'linewidth',3)
hold on
% colorvec = [0 0.6 0; 1 0.7 0.1; 0.3 0.7 0.9; 0 0 0 ; 1 0 0.1];
colorvec = [0 0.6 0; 0.3 0.7 0.9; 1 0.7 0.1; 0 0 0 ; 1 0 0.1];
for i = [1:3 5]
plot(MapBoundAve,pProbVOD_dT_SigNeg(i,:),'-o','linewidth',2,'color',...
    colorvec(i,:))
hold on
end
% plot(MapBoundAve,pProbVOD_dT_SigNeg([1:3 5],:),'--o','linewidth',2)
% ylabel([{'Significant Plant'}; {'Drying Effect (% Area)'}])
% ylabel([{'Significant Influence'}; {'on VOD (% Area)'}])
% ylabel('Areal Prevalence (%)')
% ylabel([{'Interaction'};{'Areal Prevalence (%)'}])
% ylabel([{'Percent Area with Statistically'};{'Significant Interaction (%)'}])
ylabel([{'Prevalence of Statistically'};{'Significant Interaction (% Area)'}])
legend('dT_{t-1}\rightarrow VOD_{t} (Direct)',...
    'VOD_{t-2}\rightarrow dT_{t-1}\rightarrow VOD_{t}',...
    '      \theta_{t-2}\rightarrow dT_{t-1}\rightarrow VOD_{t}',...
    '   Rs_{t-2}\rightarrow dT_{t-1}\rightarrow VOD_{t}',...
    'VPD_{t-2}\rightarrow dT_{t-1}\rightarrow VOD_{t}','fontsize',12,'position',[0.302 0.29 0.1 0.1])
% ylabel([{'Plant Drying Influence'}; {'via dT (% Area)'}])
ylim([0 70])
grid on
set(gca,'fontsize',14)
xlabel('Mean Annual Precipitation (mm)')
% title('dT \rightarrow VOD')
% text(2200,90,[{'Influence of dT Interactions on Plant Water Content'}],...
%     'fontsize',20,'horizontalalignment', 'center','FontWeight','bold')
text(2300,77,[{'Effect of dT Interactions on Plant Water Content'}],...
    'fontsize',17,'horizontalalignment', 'center','FontWeight','bold')

iplotB = subplot(2,2,2);
% hAx = axes;
% set(iplotB, 'box','off')

VOD_SM_Frac1 = reshape(VOD_SM_Frac,[size(VOD_SM_Frac,1)...
    *size(VOD_SM_Frac,2) size(VOD_SM_Frac,3)]);
nVOD_SM_Frac1 = VOD_SM_Frac1; nVOD_SM_Frac1(isfinite(nVOD_SM_Frac1))=1;
nVOD_SM_Frac1 = nansum(nVOD_SM_Frac1,1);
h1=boxplot(VOD_SM_Frac1(:,[1 3:5]), 'orientation', 'horizontal');
% h1=boxplot(VOD_SM_Frac1(:,[1 3:5]));
set(h1,{'linew'},{1},{'color'},{'k'},'linestyle','-')
% ylim([0 0.80])
xlim([-0.85 0.65])
% xlim([-0.65 0.65])
ylim([-0.5 5.5])

colorvec = flipud([0 0.6 0; 1 0.7 0.1; 0.3 0.7 0.9; 1 0 0.1]);
h=findobj(gca,'tag','Outliers');
delete(h)
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),colorvec(j,:),'FaceAlpha',0.4);
end

hold on 
% plot([0 0],[0 10],'--k','linewidth',2)
plot([0 0],[0.75 6],'--k','linewidth',1)
set(gca, 'Ydir', 'reverse')

YDist = -0.8; YRot = 0; YSize = 12;
set(gca,'YTickLabel',{''})
% text(-1,1,[{'VOD_{t-2}\rightarrow \theta_{t-1}\rightarrow'}; {'\rightarrow VOD_{t}'}]...
%     ,'horizontalalignment', 'center','fontsize',10)
% text(-1,2,[{'Rs_{t-2}\rightarrow \theta_{t-1}\rightarrow'}; {'\rightarrow VOD_{t}'}],'horizontalalignment', 'center','fontsize',10)
% text(-1,3,[{'dT_{t-2}\rightarrow \theta_{t-1}\rightarrow'}; {'\rightarrow VOD_{t}'}],'horizontalalignment', 'center','fontsize',10)
% text(-1,4,[{'VPD_{t-2}\rightarrow \theta_{t-1}\rightarrow'}; {'\rightarrow VOD_{t}'}],'horizontalalignment', 'center','fontsize',10)
text(YDist,1,[{'VOD_{t-2}\rightarrow \theta_{t-1}\rightarrow VOD_{t}'}]...
    ,'horizontalalignment', 'center','fontsize',YSize,'rotation',YRot)
text(YDist,2,[{'Rs_{t-2}\rightarrow \theta_{t-1}\rightarrow VOD_{t}'}],'horizontalalignment', 'center','fontsize',YSize,'rotation',YRot)
text(YDist,3,[{'dT_{t-2}\rightarrow \theta_{t-1}\rightarrow VOD_{t}'}],'horizontalalignment', 'center','fontsize',YSize,'rotation',YRot)
text(YDist,4,[{'VPD_{t-2}\rightarrow \theta_{t-1}\rightarrow VOD_{t}'}],'horizontalalignment', 'center','fontsize',YSize,'rotation',YRot)
% text(-0.8,1,[{'X = VOD'}],'horizontalalignment', 'center','fontsize',14)
% text(-0.8,2,[{'X = Rs'}],'horizontalalignment', 'center','fontsize',14)
% text(-0.8,3,[{'X = dT'}],'horizontalalignment', 'center','fontsize',14)
% text(-0.8,4,[{'X = VPD'}],'horizontalalignment', 'center','fontsize',14)
set(gca,'YColor','none')

% set(gca,'XTickLabel',[{'VOD\rightarrow \theta\rightarrow VOD'}; ...
%     {'Rs\rightarrow \theta\rightarrow VOD'};...
%     {'dT\rightarrow \theta\rightarrow VOD'};...
%     {'VPD\rightarrow \theta\rightarrow VOD'}],'fontsize',9)
% set(gca,'XTickLabel',['A';'A';'A';'A'],'fontsize',18)
% xlabel([{'Relative Magnitude of'}; {'Indirect Influence on VOD'}])
xlabel(['Interaction Relative Magnitude'])
% ylabel([{'Significant VOD'}; {'Influence (% Area)'}])
hAx=gca;
hAx.XAxis.TickLabelInterpreter='tex';
set(gca,'fontsize',14)
% xtickangle(20)
% text(0.5,4.5,'$\frac{{\Delta}VOD_{X{\rightarrow}{\theta}{\rightarrow}VOD}}{{\Delta}VOD_{{\theta}{\rightarrow}VOD}}$',...
%     'Interpreter','latex','fontsize',20,'horizontalalignment', 'center')
text(0,0.1,'$\frac{{\Delta}VOD_{X_{t-2}{\rightarrow}{{\theta}_{t-1}}{\rightarrow}VOD_{t}}}{{\Delta}VOD_{{{\theta}_{t-1}}{\rightarrow}VOD_{t}}}$',...
    'Interpreter','latex','fontsize',19,'horizontalalignment', 'center')
text(0.45,5,[{'Enhances'}; {'Plant Drying'}],'fontsize',13,'horizontalalignment', 'center')
text(-0.45,5,[{'Suppresses'}; {'Plant Drying'}],'fontsize',13,'horizontalalignment', 'center')
% annotation('textarrow',[0.64 0.64],[0.85 0.91])
% annotation('textarrow',[0.64 0.64],[0.83 0.77])
annotation('textarrow',[0.815 0.865],[0.94 0.94]-0.33)
annotation('textarrow',[0.76 0.71]+0.04,[0.94 0.94]-0.33)

grid on

%%%%% dT Boxplot %%%%
iplotD = subplot(2,2,4);
VOD_dT_Frac1 = reshape(VOD_dT_Frac,[size(VOD_dT_Frac,1)...
    *size(VOD_dT_Frac,2) size(VOD_dT_Frac,3)]);
nVOD_dT_Frac1 = VOD_dT_Frac1; nVOD_dT_Frac1(isfinite(nVOD_dT_Frac1))=1;
nVOD_dT_Frac1 = nansum(nVOD_dT_Frac1,1);
h1=boxplot(VOD_dT_Frac1(:,[1:3 5]), 'orientation', 'horizontal');
set(h1,{'linew'},{1},{'color'},{'k'},'linestyle','-')
xlim([-0.85 0.65])
ylim([-0.5 5.5])
h=findobj(gca,'tag','Outliers');
delete(h)
h = findobj(gca,'Tag','Box');
colorvec = flipud([0 0.6 0; 0.3 0.7 0.9; 1 0.7 0.1 ; 1 0 0.1]);
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),colorvec(j,:),'FaceAlpha',0.4);
end
hold on 
% text(1,-1,{'VOD\rightarrow dT\rightarrow VOD'},'horizontalalignment', 'center','fontsize',12,'rotation',20)
% text(2,-1,{'\theta\rightarrow dT\rightarrow VOD'},'horizontalalignment', 'center','fontsize',12,'rotation',20)
% text(3,-1,{'Rs\rightarrow dT\rightarrow VOD'},'horizontalalignment', 'center','fontsize',12,'rotation',20)
% text(4,-1,{'VPD\rightarrow dT\rightarrow VOD'},'horizontalalignment', 'center','fontsize',12,'rotation',20)
% text(-0.95,1,{'VOD_{t-2}\rightarrow dT_{t-1}\rightarrow VOD_{t}'},'horizontalalignment', 'center','fontsize',10)
% text(-0.95,2,{'\theta_{t-2}\rightarrow dT_{t-1}\rightarrow VOD_{t}'},'horizontalalignment', 'center','fontsize',10)
% text(-0.95,3,{'Rs_{t-2}\rightarrow dT_{t-1}\rightarrow VOD_{t}'},'horizontalalignment', 'center','fontsize',10)
% text(-0.95,4,{'VPD_{t-2}\rightarrow dT_{t-1}\rightarrow VOD_{t}'},'horizontalalignment', 'center','fontsize',10)
% text(-0.8,1,[{'X = VOD'}],'horizontalalignment', 'center','fontsize',14)
% text(-0.8,2,[{'X = \theta'}],'horizontalalignment', 'center','fontsize',14)
% text(-0.8,3,[{'X = Rs'}],'horizontalalignment', 'center','fontsize',14)
% text(-0.8,4,[{'X = VPD'}],'horizontalalignment', 'center','fontsize',14)

text(YDist,1 - 0.25,[{'VOD_{t-2}\rightarrow dT_{t-1}\rightarrow VOD_{t}'}]...
    ,'horizontalalignment', 'center','fontsize',YSize,'rotation',YRot)
text(YDist,2- 0.25,[{'\theta_{t-2}\rightarrow dT_{t-1}\rightarrow VOD_{t}'}],'horizontalalignment', 'center','fontsize',YSize,'rotation',YRot)
text(YDist,3- 0.25,[{'Rs_{t-2}\rightarrow dT_{t-1}\rightarrow VOD_{t}'}],'horizontalalignment', 'center','fontsize',YSize,'rotation',YRot)
text(YDist,4- 0.25,[{'VPD_{t-2}\rightarrow dT_{t-1}\rightarrow VOD_{t}'}],'horizontalalignment', 'center','fontsize',YSize,'rotation',YRot)

hold on 
set(gca, 'Ydir', 'reverse')
% plot([-1 10],[0 0],'--k','linewidth',2)
% plot([0 0],[-1 5.2],'--k','linewidth',1)
plot([0 0],[0.75 6],'--k','linewidth',1)
set(gca,'YTickLabel',{''})
xlabel(['Interaction Relative Magnitude'])
set(gca,'YColor','none')
% text(0.5,0.4,[{'Enhances'}; {'Plant Drying'}],'fontsize',12,'horizontalalignment', 'center')
% text(-0.5,0.4,[{'Suppresses'}; {'Plant Drying'}],'fontsize',12,'horizontalalignment', 'center')
% annotation('textarrow',[0.64 0.64],[0.52 0.58])
% annotation('textarrow',[0.64 0.64],[0.5 0.44])
% annotation('textarrow',[0.79 0.84],[0.61 0.61])
% annotation('textarrow',[0.77 0.72],[0.61 0.61])
% annotation('textarrow',[0.64 0.64],[0.5 0.44])
grid on
text(0.45,5,[{'Enhances'}; {'Plant Drying'}],'fontsize',13,'horizontalalignment', 'center')
text(-0.45,5,[{'Suppresses'}; {'Plant Drying'}],'fontsize',13,'horizontalalignment', 'center')
annotation('textarrow',[0.815 0.865],[0.94 0.94]-0.83)
annotation('textarrow',[0.76 0.71]+0.04,[0.94 0.94]-0.83)

hAx=gca;
hAx.XAxis.TickLabelInterpreter='tex';
set(gca,'fontsize',14)
% xtickangle(20)
text(0,0.1,'$\frac{{\Delta}VOD_{X_{t-2}{\rightarrow}{dT_{t-1}}{\rightarrow}VOD_{t}}}{{\Delta}VOD_{{dT_{t-1}}{\rightarrow}VOD_{t}}}$',...
    'Interpreter','latex','fontsize',19,'horizontalalignment', 'center')

set(gcf,'position',[10 10 650 550])

AMSA = get(iplotA,'Position'); % top right fig
set(iplotA,'Position',[AMSA(1)-0.035 AMSA(2)-0.01 AMSA(3)+0.06 AMSA(4)+0.02]);
AMSA = get(iplotA,'Position'); % top right fig

AMSB = get(iplotB,'Position'); % top right fig
set(iplotB,'Position',[AMSB(1)+0.02 AMSA(2) AMSA(3) AMSA(4)]);
AMSB = get(iplotB,'Position'); % top right fig

AMSC = get(iplotC,'Position'); % top right fig
set(iplotC,'Position',[AMSA(1) AMSC(2)-0.04 AMSA(3) AMSA(4)]);
AMSC = get(iplotC,'Position'); % top right fig

AMSD = get(iplotD,'Position'); % top right fig
set(iplotD,'Position',[AMSB(1) AMSC(2) AMSB(3) AMSA(4)]);
AMSD = get(iplotD,'Position'); % top right fig

annotation('line',[0 1],[0.5 0.5],'linewidth',1.5,'linestyle','--')
annotation('line',[0 1],[0 0],'linewidth',3)
annotation('line',[0 0],[0 1],'linewidth',3)
annotation('line',[1 1],[0 1],'linewidth',3)
annotation('line',[0 1],[1 1],'linewidth',3)

Atext = ['A'];
annotation('textbox',[0.005 0.905 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',24,'EdgeColor','none')
Atext = ['B'];
annotation('textbox',[0.005 0.405 0.1 0.1],'string',...
    Atext,'FitBoxToText','on','fontsize',24,'EdgeColor','none')
