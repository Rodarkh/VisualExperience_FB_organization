function generate_Figure_5(initial_LM_dRF, initial_V1_dRF, dendriticSize, V1_RFshape, connectivityRule, borderOfEllipseFromSD, survivalGain, nSamples, n_iter) 
% Model created and implemented by Rodrigo F. Dias and Leopoldo Petreanu
% Code written by Rodrigo F. Dias and Leopoldo Petreanu - Cortical Circuits Laboratory - Champalimaud Research - 2024
%
% Work under revision (NEURON-D-22-01431R1)
% "Visual experience reduces the spatial redundancy between cortical feedback inputs and primary visual cortex neurons." 
%
%% settings
currSettings = ['dend ' ,num2str(dendriticSize), '; surv ',num2str(survivalGain), '; ',num2str(connectivityRule), '; SD ',num2str(borderOfEllipseFromSD)];

groups     = {'DRP0'        ,'DRP21'       ,'NR','NR_L5'    ,'DRP21_L5'   ,'NR_L23'          ,'DRP21_L23'};
color2plot = {rgb('DarkRed'),rgb('Salmon') ,'k' ,rgb('Navy'),rgb('Sienna'),rgb('DeepSkyBlue'),rgb('DarkOrange')};

color2plot_init = [204 51 159]/256;

%%
nGroups=length(groups);
dRF_amp_counts_init_all=[];
vert_data_init_all=[];
meanAreaInBin_init_all=[];
dRF_amp_init_all=[];
% Loading data save when running model for each group, all iterations
for iii=1:nGroups
    LM_Group1_RFs{iii} = [groups{iii},'_OS90'];
    LM_Group2_RFs{iii} = [groups{iii},'_OS0'];

    pathLoad = [ 'Model', filesep, currSettings, filesep, groups{iii}, filesep, 'analysis'];

    clear x xVars
    for iter=1:n_iter
        x{iter}=load([pathLoad filesep 'data_',initial_LM_dRF,'_',initial_V1_dRF,'_G1_',groups{iii},'_OS90_G2_',groups{iii},'_OS0_',V1_RFshape,'_',num2str(nSamples/1000),'k__iter',num2str(iter),'.mat']);
    end
    xVars= fieldnames(x{1});

    for i=1:length(xVars)
        for iter=1:n_iter
            if strcmp(xVars{i},'survivalG1')
                survivalG1{iii}(:,:,iter) = x{iter}.(xVars{i});
            elseif strcmp(xVars{i},'survivalG2')
                survivalG2{iii}(:,:,iter) = x{iter}.(xVars{i});
            elseif strcmp(xVars{i},'counts_G1')
                counts_G1{iii}(:,:,iter) = x{iter}.(xVars{i});
            elseif strcmp(xVars{i},'counts_G2')
                counts_G2{iii}(:,:,iter) = x{iter}.(xVars{i});
            elseif strcmp(xVars{i},'counts_V1_G1')
                counts_V1_G1{iii}(:,:,iter) = x{iter}.(xVars{i});
            elseif strcmp(xVars{i},'counts_V1_G2')
                counts_V1_G2{iii}(:,:,iter) = x{iter}.(xVars{i});
            elseif strcmp(xVars{i},'dRF_amp_counts_bins')
                dRF_amp_counts_bins(:,:,iter) = x{iter}.(xVars{i});
            elseif strcmp(xVars{i},'binMeans_Group1_plot')
                binMeans_Group1_plot{iii}(:,:,iter) = x{iter}.(xVars{i});
            elseif strcmp(xVars{i},'binMeans_Group2_plot')
                binMeans_Group2_plot{iii}(:,:,iter) = x{iter}.(xVars{i});
            elseif strcmp(xVars{i},'vert_data_surv')
                vert_data_surv{iii}(:,:,iter) = x{iter}.(xVars{i});
            elseif strcmp(xVars{i},'vert_data_init')
                vert_data_init{iii}(:,:,iter) = x{iter}.(xVars{i});
            elseif strcmp(xVars{i},'dRF_amp_counts_surv')
                dRF_amp_counts_surv{iii}(:,:,iter) = x{iter}.(xVars{i});
            elseif strcmp(xVars{i},'dRF_amp_counts_init')
                dRF_amp_counts_init{iii}(:,:,iter) = x{iter}.(xVars{i});
            elseif strcmp(xVars{i},'rhoG1')
                rhoG1{iii}(:,:,iter) = x{iter}.(xVars{i});
            elseif strcmp(xVars{i},'rhoG2')
                rhoG2{iii}(:,:,iter) = x{iter}.(xVars{i});
            elseif strcmp(xVars{i},'rhoSurvivingG1')
                rhoSurvivingG1{iii,iter}(:,:) = x{iter}.(xVars{i});
            elseif strcmp(xVars{i},'rhoSurvivingG2')
                rhoSurvivingG2{iii,iter}(:,:) = x{iter}.(xVars{i});
            elseif strcmp(xVars{i},'survival_group1')
                survival_group1{iii}(:,:,iter) = x{iter}.(xVars{i});
            elseif strcmp(xVars{i},'survival_group2')
                survival_group2{iii}(:,:,iter) = x{iter}.(xVars{i});
            end
        end
    end

    % concatenation of data when needed to average across the two
    % orientations
    dRF_amp_counts_bins = squeeze(dRF_amp_counts_bins(:,:,1));

    dRF_amp_init_all_perGroup{iii}(:,1,1) = nanmedian(cat(1,rhoG1{iii},rhoG2{iii})); 
    for iter=1:n_iter
        dRF_amp_surv_all{iii}(iter,1) = nanmedian(cat(1,rhoSurvivingG1{iii,iter},rhoSurvivingG2{iii,iter}));
    end

    dRF_amp_init_all = cat(1,dRF_amp_init_all,dRF_amp_init_all_perGroup{iii});
    
    vert_data_init_all = cat(3,vert_data_init_all,vert_data_init{iii});
    dRF_amp_counts_init_all = cat(3,dRF_amp_counts_init_all,dRF_amp_counts_init{iii});
end

%%
% (radial) corrected dRF histograms
rho=1:1:dRF_amp_counts_bins(end);
for r=1:1:dRF_amp_counts_bins(end)
    area_per_rho(r) = (pi*r^2) -  (pi*(r-1)^2);
end
for iii=1:nGroups
    dRF_amp_counts_surv_corr{iii} = dRF_amp_counts_surv{iii}./area_per_rho;
    dRF_amp_counts_init_corr{iii} = dRF_amp_counts_init{iii}./area_per_rho;
end
dRF_amp_counts_init_all_corr = dRF_amp_counts_init_all./area_per_rho;

customText=[];
saveFigFolder = ['Model', filesep, currSettings];
saveFinalFigFolder = ['Model'];
titleText = [currSettings, '; ',num2str(nSamples/1000),'k samples'];

% Figure display settings
mainColorMap = 'hot';
divColorMap = cbrewer('div','RdBu',31);
min_counts = 30;

OverlapFractionClims = [-.075 .075];
Maps2DLimit = 47.5;
Maps2DTick = -60:30:60;

%% Data per Group 
% Average of results across all iterations of the model for each group
for iii=1:nGroups

    figure
    xPlots = 3;

    % initial dRF
    ax1=subplot(2,xPlots,1)
    h=imagesc((squeeze(nanmean(counts_G1{iii},3))'));
    hold on;
    dx=2.5; 
    XD=-80+dx:dx:80-dx;
    YD=-80+dx:dx:80-dx;
    hold on; plot(0,0,'k+','MarkerSize',5);
    for jj=[10] %Draw stimuli circles
        th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
        plot(xunit, yunit,'--k','LineWidth',1);
    end
    set(h,'Xdata',XD,'Ydata',YD );
    c=colorbar;
    set(gca,'YDir','normal','xtick',Maps2DTick,'xlim',[-Maps2DLimit, Maps2DLimit],'ytick',Maps2DTick,'ylim',[-Maps2DLimit, Maps2DLimit])
    title(['initial \DeltaRF dist. from ',initial_LM_dRF, ' for ', LM_Group1_RFs{iii}]);
    axis square

    ax4=subplot(2,xPlots,4)
    h=imagesc((squeeze(nanmean(counts_G2{iii},3))'));
    hold on;
    dx=2.5;
    XD=-80+dx:dx:80-dx;
    YD=-80+dx:dx:80-dx;
    hold on; plot(0,0,'k+','MarkerSize',5);
    for jj=[10] %Draw stimuli circles
        th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
        plot(xunit, yunit,'--k','LineWidth',1);
    end
    set(h,'Xdata',XD,'Ydata',YD );
    c=colorbar;
    set(gca,'YDir','normal','xtick',Maps2DTick,'xlim',[-Maps2DLimit, Maps2DLimit],'ytick',Maps2DTick,'ylim',[-Maps2DLimit, Maps2DLimit])
    title(['initial \DeltaRF dist. from ',initial_LM_dRF, ' for ', LM_Group2_RFs{iii}]);
    axis square

    % survival dRF Group 1
    ax2=subplot(2,xPlots,2)
    h=imagesc(nanmean(survivalG1{iii},3));
    hold on;
    dx=2.5; 
    XD=-80+dx:dx:80-dx;
    YD=-80+dx:dx:80-dx;
    hold on; plot(0,0,'k+','MarkerSize',5);
    for jj=[10] %Draw stimuli circles
        th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
        plot(xunit, yunit,'--k','LineWidth',1);
    end
    set(h,'Xdata',XD,'Ydata',YD );
    c=colorbar;
    set(gca,'YDir','normal','xtick',Maps2DTick,'xlim',[-Maps2DLimit, Maps2DLimit],'ytick',Maps2DTick,'ylim',[-Maps2DLimit, Maps2DLimit])
    title(['Surviving \DeltaRF dist. - ' LM_Group1_RFs{iii}  ]);
    axis square

    % survival dRF Group 2
    ax5=subplot(2,xPlots,5)
    h=imagesc(nanmean(survivalG2{iii},3));
    hold on;
    dx=2.5; 
    XD=-80+dx:dx:80-dx;
    YD=-80+dx:dx:80-dx;
    hold on; plot(0,0,'k+','MarkerSize',5);
    for jj=[10] %Draw stimuli circles
        th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
        plot(xunit, yunit,'--k','LineWidth',1);
    end
    set(h,'Xdata',XD,'Ydata',YD );
    c=colorbar;
    set(gca,'YDir','normal','xtick',Maps2DTick,'xlim',[-Maps2DLimit, Maps2DLimit],'ytick',Maps2DTick,'ylim',[-Maps2DLimit, Maps2DLimit])
    hold off
    title(['Surviving \DeltaRF dist. - ' LM_Group2_RFs{iii}  ]);
    axis square


    ax3=subplot(2,xPlots,3)
    h=imagesc(nanmean([binMeans_Group1_plot{iii}-binMeans_Group2_plot{iii}],3),OverlapFractionClims);
    hold on;
    dx=2.5; 
    XD=-80+dx:dx:80-dx;
    YD=-80+dx:dx:80-dx;
    hold on; plot(0,0,'k+','MarkerSize',5);
    for jj=[10] %Draw stimuli circles
        th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
        plot(xunit, yunit,'--k','LineWidth',1);
    end
    set(h,'Xdata',XD,'Ydata',YD );
    c=colorbar;
    set(gca,'YDir','normal','xtick',Maps2DTick,'xlim',[-Maps2DLimit, Maps2DLimit],'ytick',Maps2DTick,'ylim',[-Maps2DLimit, Maps2DLimit])
    c.Label.String = ['\DeltaOverlap counts>', num2str(min_counts)];
    title({['Mean Fraction of RF Overlap difference in \DeltaRF space'] ;[LM_Group1_RFs{iii} ' - ' LM_Group2_RFs{iii}  ]});
    axis square

    ax6=subplot(2,xPlots,6)
    currMaxValue=max(abs(nanmean(survivalG1{iii},3)/max(nanmean(survivalG1{iii},3),[],'all')-nanmean(survivalG2{iii},3)/max(nanmean(survivalG2{iii},3),[],'all')),[],'all');
    h=imagesc(nanmean(survivalG1{iii},3)/max(nanmean(survivalG1{iii},3),[],'all')-nanmean(survivalG2{iii},3)/max(nanmean(survivalG2{iii},3),[],'all'),[-currMaxValue currMaxValue]);
    hold on;
    plot(0,0,'k+','MarkerSize',5);
    dx=2.5; 
    XD=-80+dx:dx:80-dx;
    YD=-80+dx:dx:80-dx;
    hold on; plot(0,0,'k+','MarkerSize',5);
    for jj=[10] %Draw stimuli circles
        th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
        plot(xunit, yunit,'--k','LineWidth',1);
    end
    set(h,'Xdata',XD,'Ydata',YD );
    c=colorbar;
    set(gca,'YDir','normal','xtick',Maps2DTick,'xlim',[-Maps2DLimit, Maps2DLimit],'ytick',Maps2DTick,'ylim',[-Maps2DLimit, Maps2DLimit])
    c.Label.String = '\DeltaDensity';
    title({ ['Surviving \DeltaRF distribution'];[LM_Group1_RFs{iii} ' - ' LM_Group2_RFs{iii}  ]});
    axis square

    colormap(ax1,mainColorMap)
    colormap(ax2,mainColorMap)
    colormap(ax4,mainColorMap)
    colormap(ax5,mainColorMap)
    colormap(ax3,divColorMap)
    colormap(ax6,divColorMap)

    sgtitle({[groups{iii},': Average across 10 iterations. ' ,titleText]; ['\DeltaRF starting distribution: ',initial_LM_dRF,'; ',LM_Group1_RFs{iii},' vs ',LM_Group2_RFs{iii},'; ', connectivityRule, ' connectivity rule to V1; survival rate exponential factor: ',num2str(survivalGain), '; fraction of LMV1 joint RF']},'color',color2plot{iii})
    set(gcf,'Position', [10         10        1940         920]);

    saveas(gcf,[saveFigFolder filesep 'Simulation_Reps_Average_',groups{iii}],'fig')
    saveas(gcf,[saveFigFolder filesep 'Simulation_Reps_Average_',groups{iii}],'svg')
    export_fig([saveFigFolder filesep 'Simulation_Reps_Average_',groups{iii}], '-png')
end

close all
%% Figure 5
xPlots=12;
yPlots=3;

currIter=1; % for single example plots C and E
iii=4; %NR L5 

figure(5)
clear ax1 ax2 ax3 ax4 ax5 ax6

% Pannel C, top
% V1 initial
ax1=subplot(yPlots,xPlots, xPlots*0 + [1 2 ] )
h=imagesc((squeeze(counts_V1_G2{iii}(:,:,currIter)./ 5000000)'));
dx=2.5;
XD=-80+dx:dx:80-dx;
YD=-80+dx:dx:80-dx;
hold on; plot(0,0,'k+','MarkerSize',5);
for jj=[10] %Draw stimuli circles
    th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
    plot(xunit, yunit,'--w','LineWidth',1);
end
set(h,'Xdata',XD,'Ydata',YD );
c=colorbar;
set(gca,'YDir','normal','xtick',Maps2DTick,'xlim',[-Maps2DLimit, Maps2DLimit],'ytick',Maps2DTick,'ylim',[-Maps2DLimit, Maps2DLimit])
title(['initial V1 \DeltaRF dist.']);
axis square

% Pannel E, top
% survival dRF Group 1
ax2=subplot(yPlots,xPlots, xPlots*0 + [7 8] )
h=imagesc(squeeze(survivalG1{iii}(:,:,currIter)./sum(survivalG1{iii}(:,:,currIter),'all')));
dx=2.5; 
XD=-80+dx:dx:80-dx;
YD=-80+dx:dx:80-dx;
hold on; plot(0,0,'k+','MarkerSize',5);
for jj=[10] %Draw stimuli circles
    th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
    plot(xunit, yunit,'--k','LineWidth',1);
end
set(h,'Xdata',XD,'Ydata',YD );
c=colorbar;
c.Label.String = 'Fraction';
set(gca,'YDir','normal','xtick',Maps2DTick,'xlim',[-Maps2DLimit, Maps2DLimit],'ytick',Maps2DTick,'ylim',[-Maps2DLimit, Maps2DLimit])
title(['Surviving LM \DeltaRF - ' LM_Group1_RFs{iii}  ]);
axis square

% Pannel C, bottom
% initial LM dRF for Group2
ax3=subplot(yPlots,xPlots, xPlots*1 + [1 2 ] )
h=imagesc((squeeze(counts_G2{iii}(:,:,currIter)./sum(counts_G2{iii}(:,:,currIter),'all'))'));
dx=2.5; 
XD=-80+dx:dx:80-dx;
YD=-80+dx:dx:80-dx;
hold on; plot(0,0,'k+','MarkerSize',5);
for jj=[10] %Draw stimuli circles
    th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
    plot(xunit, yunit,'--k','LineWidth',1);
end
set(h,'Xdata',XD,'Ydata',YD );
c=colorbar;
set(gca,'YDir','normal','xtick',Maps2DTick,'xlim',[-Maps2DLimit, Maps2DLimit],'ytick',Maps2DTick,'ylim',[-Maps2DLimit, Maps2DLimit])
title(['initial LM \DeltaRF dist. from ',initial_LM_dRF, ' for ', LM_Group2_RFs{iii}]);
axis square

% Pannel E, bottom
% survival LM dRF for Group 2
ax4=subplot(yPlots,xPlots, xPlots*1 + [7 8] )
h=imagesc((squeeze(survivalG2{iii}(:,:,currIter)./sum(survivalG2{iii}(:,:,currIter),'all'))));    
dx=2.5;
XD=-80+dx:dx:80-dx;
YD=-80+dx:dx:80-dx;
hold on; plot(0,0,'k+','MarkerSize',5);
for jj=[10] %Draw stimuli circles
    th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
    plot(xunit, yunit,'--k','LineWidth',1);
end
set(h,'Xdata',XD,'Ydata',YD );
c=colorbar;
c.Label.String = 'Fraction';

set(gca,'YDir','normal','xtick',Maps2DTick,'xlim',[-Maps2DLimit, Maps2DLimit],'ytick',Maps2DTick,'ylim',[-Maps2DLimit, Maps2DLimit])
title(['Surviving LM \DeltaRF - ' LM_Group2_RFs{iii}  ]);
axis square

% Pannel D
% Difference of Fraction of RF overlap per dRF bin across the two initial
% distributions (before applying survival rule)
ax5=subplot(yPlots,xPlots, [xPlots*0+4,xPlots*0+5, xPlots*0+6 , xPlots*1+4, xPlots*1+5, xPlots*1+6] )
h=imagesc(nanmean([binMeans_Group1_plot{iii}-binMeans_Group2_plot{iii}],3),OverlapFractionClims);
dx=2.5; 
XD=-80+dx:dx:80-dx;
YD=-80+dx:dx:80-dx;
hold on; plot(0,0,'k+','MarkerSize',5);
for jj=[10] %Draw stimuli circles
    th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
    plot(xunit, yunit,'--k','LineWidth',1);
end
set(h,'Xdata',XD,'Ydata',YD );
c=colorbar;
set(gca,'YDir','normal','xtick',Maps2DTick,'xlim',[-Maps2DLimit, Maps2DLimit],'ytick',Maps2DTick,'ylim',[-Maps2DLimit, Maps2DLimit])
c=colorbar;
title({['\DeltaOverlap (',LM_Group1_RFs{iii} ' - ' LM_Group2_RFs{iii} ,')']});
axis square

% Pannel F
% Difference of fraction of inputs per dRF bin across the two surviving
% distributions
ax6=subplot(yPlots,xPlots, [xPlots*0+10,xPlots*0+11, xPlots*0+12 , xPlots*1+10, xPlots*1+11, xPlots*1+12] )
currMaxValue=max(abs(nanmean(survivalG1{iii},3)/max(nanmean(survivalG1{iii},3),[],'all')-nanmean(survivalG2{iii},3)/max(nanmean(survivalG2{iii},3),[],'all')),[],'all');
h=imagesc(nanmean(survivalG1{iii},3)/max(nanmean(survivalG1{iii},3),[],'all')-nanmean(survivalG2{iii},3)/max(nanmean(survivalG2{iii},3),[],'all'),[-currMaxValue currMaxValue]);
dx=2.5; %160/64
XD=-80+dx:dx:80-dx;
YD=-80+dx:dx:80-dx;
hold on; plot(0,0,'k+','MarkerSize',5);
for jj=[10] %Draw stimuli circles
    th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
    plot(xunit, yunit,'--k','LineWidth',1);
end
set(h,'Xdata',XD,'Ydata',YD );
c=colorbar;
set(gca,'YDir','normal','xtick',Maps2DTick,'xlim',[-Maps2DLimit, Maps2DLimit],'ytick',Maps2DTick,'ylim',[-Maps2DLimit, Maps2DLimit])
c=colorbar;
title({ ['\Delta Fraction of inputs (',LM_Group1_RFs{iii} ' - ' LM_Group2_RFs{iii}  ,')']});
axis square

% Pannel G, left
% Comparison of initial dRF distribution vs surviving
subplot(yPlots,xPlots,xPlots*2 + [1 ])
clear h0 h1
hold on
iii=4; % NR L5
plotLWEA(dRF_amp_counts_bins(1:end-1),nanmean(dRF_amp_counts_surv{iii}(1,:,:),3),nanstd(dRF_amp_counts_surv{iii}(1,:,:),0,3),color2plot{iii},0.3,[])
plotLWEA(dRF_amp_counts_bins(1:end-1),nanmean(dRF_amp_counts_init{iii}(1,:,:),3),nanstd(dRF_amp_counts_init{iii}(1,:,:),0,3),color2plot_init,0.3,[])
h0=plot(dRF_amp_counts_bins(1:end-1),nanmean(dRF_amp_counts_init{iii}(1,:,:),3),'-','Color',color2plot_init,'LineWidth',2);
h1=plot(dRF_amp_counts_bins(1:end-1),nanmean(dRF_amp_counts_surv{iii}(1,:,:),3),'-','Color',color2plot{iii},'LineWidth',2);
xlim([dRF_amp_counts_bins(1) dRF_amp_counts_bins(end)])
ylabel('Fraction')
xlabel('||\DeltaRF||(º)')
set(gca,'xtick',0:10:100,'xlim',[0 41],'ytick',0:0.02:1);
hold off

% Pannel G, right
% Comparison of initial dRF distribution (radial) vs surviving
subplot(yPlots,xPlots,xPlots*2 + [2 ])
hold on
iii=4; % NR L5
plotLWEA(dRF_amp_counts_bins(1:end-1),nanmean(dRF_amp_counts_surv_corr{iii}(1,:,:),3),nanstd(dRF_amp_counts_surv_corr{iii}(1,:,:),0,3),color2plot{iii},0.3,[])
plotLWEA(dRF_amp_counts_bins(1:end-1),nanmean(dRF_amp_counts_init_corr{iii}(1,:,:),3),nanstd(dRF_amp_counts_init_corr{iii}(1,:,:),0,3),color2plot_init,0.3,[])
h0=plot(dRF_amp_counts_bins(1:end-1),nanmean(dRF_amp_counts_init_corr{iii}(1,:,:),3),'-','Color',color2plot_init,'LineWidth',2);
h1=plot(dRF_amp_counts_bins(1:end-1),nanmean(dRF_amp_counts_surv_corr{iii}(1,:,:),3),'-','Color',color2plot{iii},'LineWidth',3);
xlim([dRF_amp_counts_bins(1) dRF_amp_counts_bins(end)])
ylabel('Fraction per deg^2 x 1000')
xlabel('||\DeltaRF||(º)')
set(gca,'xtick',0:10:100,'xlim',[0 41]);
hold off

% Pannel H
% mean of Median dRF per iteration comparison across groups
subplot(yPlots,xPlots,xPlots*2 + [4 5 6])
iii=0;
hold on;
plot(iii,nanmean(dRF_amp_init_all,1),'o','Color',color2plot_init,'LineWidth',2,'markerSize',8,'MarkerFaceColor',color2plot_init,'MarkerEdgeColor','none')
errorbar(iii,nanmean(dRF_amp_init_all,1),nanstd(dRF_amp_init_all,1),'Color',color2plot_init,'LineWidth',2)

for iii=1:nGroups
    plot(iii,nanmean(dRF_amp_surv_all{iii},1),'o','Color',color2plot{iii},'LineWidth',2,'markerSize',8,'MarkerFaceColor',color2plot{iii},'MarkerEdgeColor','none')
    errorbar(iii,nanmean(dRF_amp_surv_all{iii},1),nanstd(dRF_amp_surv_all{iii},1),'Color',color2plot{iii},'LineWidth',2)
end
plot([-1 10],[1 1],'k--')
set(gca,'xlim',[-.5 nGroups+.5],'xtick',0:nGroups,'xticklabel',{'init',groups{:}},'ylim',[4 20],'ytick',0:5:20)
ylabel('median ||\DeltaRF||(º)')
title('Median |\DeltaRF|')

% Pannel I
% Comparison of initial fraction of boutons deviating in Elevation vs
% surviving one
subplot(yPlots,xPlots,xPlots*2 + [ 8 9])
clear h1 h2
hold on
iii=4;
plotLWEA(1:2,nanmean(vert_data_surv{iii},3),nanstd(vert_data_surv{iii},0,3),color2plot{iii},0.3,[])
plotLWEA(1:2,squeeze(nanmean(vert_data_init{iii},3)),squeeze(nanstd(vert_data_init{iii},0,3)),color2plot_init,0.3,[])
plot(1:2,squeeze(nanmean(vert_data_init{iii},3)),'Color',color2plot_init,'LineWidth',2)
plot(1:2,squeeze(nanmean(vert_data_surv{iii},3)),'Color',color2plot{iii},'LineWidth',2)
ylabel('Fraction of boutons deviating in elevation')
set(gca,'xlim',[0.9 2.1],'xtick',1:2,'xticklabel',{'Horizontal','Vertical'},'ylim',[0.14 .38],'ytick',0:0.05:1)
hold off

% Pannel J
% mean colinear depletion index across iterations
subplot(yPlots,xPlots,xPlots*2 + [10 11 12])
hold on
iii=0;
plot(iii,nanmean(vert_data_init_all(:,1,:)-vert_data_init_all(:,2,:),3),'o','Color',color2plot_init,'LineWidth',2,'markerSize',8,'MarkerFaceColor',color2plot_init,'MarkerEdgeColor','none')
errorbar(iii,nanmean(vert_data_init_all(:,1,:)-vert_data_init_all(:,2,:),3),nanstd(vert_data_init_all(:,1,:)-vert_data_init_all(:,2,:),0,3),'Color',color2plot_init,'LineWidth',2)

for iii=1:nGroups
    plot(iii,nanmean(vert_data_surv{iii}(:,1,:)-vert_data_surv{iii}(:,2,:),3),'o','Color',color2plot{iii},'LineWidth',2,'markerSize',8,'MarkerFaceColor',color2plot{iii},'MarkerEdgeColor','none')
    errorbar(iii,nanmean(vert_data_surv{iii}(:,1,:)-vert_data_surv{iii}(:,2,:),3),nanstd(vert_data_surv{iii}(:,1,:)-vert_data_surv{iii}(:,2,:),0,3),'Color',color2plot{iii},'LineWidth',2)
end
plot([-1 10],[0 0],'k--')
set(gca,'xlim',[-.5 nGroups+.5],'xtick',0:nGroups,'xticklabel',{'init',groups{:}},'ylim',[-.0075 .075],'ytick',0:0.01:2)
title('Collinear depletion index')

set(findobj(gcf,'type','axes'),'FontSize',10)

colormap(ax1,mainColorMap)
colormap(ax2,mainColorMap)
colormap(ax3,mainColorMap)
colormap(ax4,mainColorMap)
colormap(ax5,divColorMap)
colormap(ax6,divColorMap)

set(5,'Position', [0     40        1920         940]);

saveas(gcf,[saveFinalFigFolder filesep 'Figure_5_CJ'],'fig')
saveas(gcf,[saveFinalFigFolder filesep 'Figure_5_CJ'],'svg')
export_fig([saveFinalFigFolder filesep 'Figure_5_CJ'], '-png')

end