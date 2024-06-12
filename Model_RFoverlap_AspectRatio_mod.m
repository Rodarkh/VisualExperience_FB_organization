function Model_RFoverlap_AspectRatio_mod(dRF_LM, dRF_V1,dentriticSize, LM_Group1_RFs, LM_Group2_RFs, V1_RFs, connectivityRule, borderOfEllipseFromSD, survivalGain, numberOfSamples, ecc_factor, iterationCycleID)
%% HELP:
% Model_RFoverlap_AspectRatio_mod('DRP0', 'V1', 15, 'NR_L5_OS90', 'NR_L5_OS0', 'V1', 'isoAngle',2.5,15,5000000,1.2,iteration,'SaveLocation')

% Model created and implemented by Rodrigo F. Dias and Leopoldo Petreanu
% Code written by Rodrigo F. Dias and Leopoldo Petreanu - Cortical Circuits Laboratory - Champalimaud Research - 2024
%
% Work under revision (NEURON-D-22-01431R1)
% "Visual experience reduces the spatial redundancy between cortical feedback inputs and primary visual cortex neurons." 
%
%% Sample From 2D Gaussian of dRFs
% Load parameters from experimental values. Generate dRF distributions of
% RF centers for different groups
[AziSD, EleSD, twoDimGaussianAngle]=getdRFGaussianParameters(dRF_LM);
rotated_samples_LM_Group1=sampleFromRotated2DGaussian_LP([0 0],AziSD,EleSD,-twoDimGaussianAngle,numberOfSamples);% adding MINUS before angle to have proper direction (clockwise)
rotated_samples_LM_Group2=sampleFromRotated2DGaussian_LP([0 0],AziSD,EleSD,-twoDimGaussianAngle,numberOfSamples);
customShortText = ['_iter',num2str(iterationCycleID)];

% V1 neurons dRFGaussian
[AziSD, EleSD, twoDimGaussianAngle]=getdRFGaussianParameters(dRF_V1);
rotated_samples_V1_1=sampleFromRotated2DGaussian_LP([0 0],AziSD,EleSD,-twoDimGaussianAngle,numberOfSamples);
% one set of V1 neurons centers is generated to compare with each RF group
rotated_samples_V1_2=sampleFromRotated2DGaussian_LP([0 0],AziSD,EleSD,-twoDimGaussianAngle,numberOfSamples);

%% load parameters for RF shape LM Boutons Group 1, Group 2 and V1 neurons

[RF_longAxis_RF1, RF_longAxisSD_RF1, RF_shortAxis_RF1, RF_shortAxisSD_RF1, RF_angle_RF1, RF_angle_RF1_SD]=getRFParameters(LM_Group1_RFs);

[RF_longAxis_RF2, RF_longAxisSD_RF2, RF_shortAxis_RF2, RF_shortAxisSD_RF2, RF_angle_RF2, RF_angle_RF2_SD]=getRFParameters(LM_Group2_RFs);

[RF_longAxis_V1, RF_longAxisSD_V1, RF_shortAxis_V1, RF_shortAxisSD_V1, RF_angle_V1, RF_angle_V1_SD]=getRFParameters(V1_RFs);


overlapWithV1_group1=zeros(1,numberOfSamples);
overlapWithV1_group2=zeros(1,numberOfSamples);
survival_group1=zeros(1,numberOfSamples);
survival_group2=zeros(1,numberOfSamples);
ellipseArea_group1=zeros(1,numberOfSamples);
ellipseArea_group2=zeros(1,numberOfSamples);

ecc_G1=zeros(numberOfSamples,1);
ecc_G2=zeros(numberOfSamples,1);
dRFdata_G1 = zeros(numberOfSamples,2);
dRFdata_G2 = zeros(numberOfSamples,2);

HalfVisualSpan=150;

parfor i=1:numberOfSamples

    aG1=0;
    bG1=0;
    aG1_c=0;
    bG1_c=0;
    %ellipse LM Group 1
    %repeats axis generation until long is bigger than short
    while aG1<=bG1
        aG1=0;
        bG1=0;
        while aG1<2||aG1>50  % remove extreme and negative values of RF axis
            aG1=normrnd(RF_longAxis_RF1,RF_longAxisSD_RF1); %long axis
        end
        while bG1<2||bG1>50
            bG1=normrnd(RF_shortAxis_RF1,RF_shortAxisSD_RF1);% short axis
        end

        % create corrected axes that are elongated or shortened by the ecc
        % factor. this is done in a way as to keep the RF area unchanged
        aG1_c = aG1 * ecc_factor;
        bG1_c = bG1 / ecc_factor;

        % if this change makes the long axis shorter than the short axis,
        % then make them equal.
        if aG1_c<bG1_c
            aG1_c=bG1_c+0.00000001;
        end
    end
    aG1=aG1_c;
    bG1=bG1_c;

    angleFromHorzG1=-RF_angle_RF1;  % Changing RF angle to negative to have proper convention (positive is tilted to negative elevation on positive azimuth (so 45º is aligned with a 45ºgrating)
    deg1_G1=normrnd(angleFromHorzG1,RF_angle_RF1_SD); 

    % dRF center for sampled RFs if Group1
    x_G1=rotated_samples_LM_Group1(i,1)+HalfVisualSpan;
    y_G1=rotated_samples_LM_Group1(i,2)+HalfVisualSpan;

    aG2=0;
    bG2=0;
    aG2_c=0;
    bG2_c=0;
    %ellipse LM Group 2
    while aG2<=bG2
        aG2=0;
        bG2=0;
        while aG2<2||aG2>50
            aG2=normrnd(RF_longAxis_RF2,RF_longAxisSD_RF2); %long axis
        end

        while bG2<2||bG2>50
            bG2=normrnd(RF_shortAxis_RF2,RF_shortAxisSD_RF2);% short axis
        end
        aG2_c = aG2 * ecc_factor;
        bG2_c = bG2 / ecc_factor;

        % if this change makes the long axis shorter than the short axis,
        % then make them equal.
        if aG2_c<bG2_c
            aG2_c=bG2_c+0.00000001;
        end
    end
    aG2=aG2_c;
    bG2=bG2_c;

    angleFromHorzG2=-RF_angle_RF2;  % Changing RF angle to negative to have proper convention (positive is tilted to negative elevation on positive azimuth (so 45º is aligned with a 45ºgrating)
    deg1_G2=normrnd(angleFromHorzG2,RF_angle_RF2_SD); 

    % dRF center for sampled RFs if Group2
    x_G2=rotated_samples_LM_Group2(i,1)+HalfVisualSpan;
    y_G2=rotated_samples_LM_Group2(i,2)+HalfVisualSpan;


    %ellipse V1 Neurons:
    %normal
    aV1_1=0;
    bV1_1=0;
    while aV1_1<=bV1_1 % Making sure long axis is bigger than short axis
        aV1_1=0;
        bV1_1=0;
        while aV1_1<2||aV1_1>50
            aV1_1=normrnd(RF_longAxis_V1,RF_longAxisSD_V1);
        end

        while bV1_1<2||bV1_1>50
            bV1_1=normrnd(RF_shortAxis_V1,RF_shortAxisSD_V1);
        end
    end

    aV1_2=0;
    bV1_2=0;
    while aV1_2<=bV1_2
        aV1_2=0;
        bV1_2=0;
        while aV1_2<2||aV1_2>50
            aV1_2=normrnd(RF_longAxis_V1,RF_longAxisSD_V1);
        end

        while  bV1_2<2||bV1_2>50
            bV1_2=normrnd(RF_shortAxis_V1,RF_shortAxisSD_V1);
        end
    end


    % pick V1 RF angle depending on connectivity rule
    switch connectivityRule
        case 'random'
            deg_V1_1=normrnd(RF_angle_V1,RF_angle_V1_SD);
            deg_V1_2=normrnd(RF_angle_V1,RF_angle_V1_SD);

        case 'isoAngle'
            angularJitter=5;
            deg_V1_1=deg1_G1+((rand(1)-0.5)*angularJitter*2);
            deg_V1_2=deg1_G2+((rand(1)-0.5)*angularJitter*2);

        case 'crossAngle'
            angularJitter=1;
            deg_V1_1=deg1_G1+(rand(1)*angularJitter-angularJitter/2)+90;
            deg_V1_2=deg1_G2+(rand(1)*angularJitter-angularJitter/2)+90;
    end
    %% DENDRITES

    % we randomly shifted the V1 RF center
    %  position up to 15° in azimuth and 5° in elevation, which corresponds to dendritic trees
    % with a span of up to 300 µm given the cortical gain of V1 (if dentriticSize=15)
    randTheta=(rand(1)-1)*2*pi;
    randRho=(rand(1)-1)*dentriticSize;

    [denx, deny]=pol2cart(randTheta, randRho);

    x_V1_1=rotated_samples_V1_1(i,1)+HalfVisualSpan+denx;
    y_V1_1=rotated_samples_V1_1(i,2)+HalfVisualSpan+deny/3;

    randTheta=(rand(1)-1)*2*pi;
    randRho=(rand(1)-1)*dentriticSize;

    [denx, deny]=pol2cart(randTheta, randRho);

    x_V1_2=rotated_samples_V1_2(i,1)+HalfVisualSpan+denx;
    y_V1_2=rotated_samples_V1_2(i,2)+HalfVisualSpan+deny/3;
    distance=0;

    overlapWithV1_group1(i) = calculateEllipseOverlap_LPRD(aG1, bG1, deg1_G1, x_G1, y_G1, aV1_1, bV1_1, deg_V1_1, x_V1_1, y_V1_1, distance,borderOfEllipseFromSD,HalfVisualSpan);
    overlapWithV1_group2(i) = calculateEllipseOverlap_LPRD(aG2, bG2, deg1_G2, x_G2, y_G2, aV1_2, bV1_2, deg_V1_2, x_V1_2, y_V1_2, distance,borderOfEllipseFromSD,HalfVisualSpan);


    %% Survival Rule
    survivalThreshold1=exp(-survivalGain*overlapWithV1_group1(i));
    survivalThreshold2=exp(-survivalGain*overlapWithV1_group2(i));
    survival_group1(i) = rand(1)<(survivalThreshold1);
    survival_group2(i) = rand(1)<(survivalThreshold2);

    ecc_G1(i,1) = sqrt(1 - ( (bG1.*bG1)./(aG1.*aG1) ) );
    ecc_G2(i,1) = sqrt(1 - ( (bG2.*bG2)./(aG2.*aG2) ) );

    dRFdata_G1(i,:) = [ x_G1 - x_V1_1 , y_G1 - y_V1_1];
    dRFdata_G2(i,:) = [ x_G2 - x_V1_2 , y_G2 - y_V1_2];
end
%%

% Define bin edges
numberOfBins=64;
x_edges = linspace(-80, 80, numberOfBins);
y_edges = linspace(-80, 80, numberOfBins);

min_counts = 30;

% Get counts from V1 initial distribution
[counts_V1_G1, ~, ~, binx_V1_G1, biny_V1_G1] = histcounts2(rotated_samples_V1_1(:,1), rotated_samples_V1_1(:,2),x_edges, y_edges);
[counts_V1_G2, ~, ~, binx_V1_G2, biny_V1_G2] = histcounts2(rotated_samples_V1_2(:,1), rotated_samples_V1_2(:,2),x_edges, y_edges);

% calculate Mean overlap in each bin
[counts_G1, ~, ~, binx_G1, biny_G1] = histcounts2(rotated_samples_LM_Group1(:,1), rotated_samples_LM_Group1(:,2),x_edges, y_edges);
validXY_G1 = find(counts_G1>min_counts);

[counts_G2, ~, ~, binx_G2, biny_G2] = histcounts2(rotated_samples_LM_Group2(:,1), rotated_samples_LM_Group2(:,2),x_edges, y_edges);
validXY_G2 = find(counts_G2>min_counts);

% Calculate mean in each bin
bin_means_G1 = nan(size(counts_G1));
survivalHist_G1=zeros(size(counts_G1));
bin_indices1 = false(size(binx_G1));
bin_indices1_smooth = false(size(binx_G1));

for i = 1:max(binx_G1)
    for j = 1:max(biny_G1)
        bin_indices1 = binx_G1 == i & biny_G1 == j ;
        bin_indices1_smooth = binx_G1 == i & biny_G1 == j & sum(sub2ind([size(counts_G1,1) size(counts_G1,2)],i,j)==validXY_G1)>0 & sum(sub2ind([size(counts_G2,1) size(counts_G2,2)],i,j)==validXY_G2)>0;
        bin_means_G1(i,j) = nanmean(overlapWithV1_group1(bin_indices1_smooth));
        survivalHist_G1(i,j) = sum(survival_group1(bin_indices1));         
    end
end


bin_means_G2 = nan(size(counts_G2));
survivalHist_G2=zeros(size(counts_G2));
bin_indices2 = false(size(binx_G2));
bin_indices2_smooth = false(size(binx_G2));

for i = 1:max(binx_G2)
    for j = 1:max(biny_G2)
        bin_indices2 = binx_G2 == i & biny_G2 == j;
        bin_indices2_smooth = binx_G2 == i & biny_G2 == j & sum(sub2ind([size(counts_G2,1) size(counts_G2,2)],i,j)==validXY_G2)>0 & sum(sub2ind([size(counts_G1,1) size(counts_G1,2)],i,j)==validXY_G1)>0;
        bin_means_G2(i,j) = nanmean(overlapWithV1_group2(bin_indices2_smooth));
        survivalHist_G2(i,j) = sum(survival_group2(bin_indices2)); 
    end
end

binMeans_Group1=bin_means_G1';
binMeans_Group2=bin_means_G2';
survivalG1=survivalHist_G1';
survivalG2=survivalHist_G2';


%% Plot
currSettings = ['dend ' ,num2str(dentriticSize), '; surv ',num2str(survivalGain), '; ',num2str(connectivityRule), '; SD ',num2str(borderOfEllipseFromSD)];
saveFigFolder = [ 'Model', filesep, currSettings, '_AspectRatio'];
titleText = [currSettings, '; ',num2str(numberOfSamples/1000),'k samples; ecc factor ',num2str(ecc_factor,2)];

if ~exist(saveFigFolder)
    mkdir(saveFigFolder);
end

close all
mainColorMap = 'hot';
divColorMap = cbrewer('div','RdBu',31);
init_color = [204 51 159]/256;
%% Figure with model results per iteration


figure
xPlots = 6;
% initial 2D histogram dRF for V1 for Group 1
ax1=subplot(2,xPlots,1)
h=imagesc((counts_V1_G1'));
set(gca,'YDir','normal')
dx=2.5; 
XD=-80+dx:dx:80-dx;
YD=-80+dx:dx:80-dx;
hold on; plot(0,0,'k+','MarkerSize',5);
for jj=[10] %Draw circles
    th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
    plot(xunit, yunit,'--k','LineWidth',1);
end
set(h,'Xdata',XD,'Ydata',YD );
xlim([-60+dx, 60-dx]);
ylim([-60+dx, 60-dx]);
c=colorbar;
title(['initial \DeltaRF dist. from ',dRF_V1, ' for ', LM_Group1_RFs],'FontSize',8);
axis square

% initial 2D histogram dRF for V1 for Group 2
ax7=subplot(2,xPlots,xPlots+1)
h=imagesc((counts_V1_G2'));
set(gca,'YDir','normal')
dx=2.5; 
XD=-80+dx:dx:80-dx;
YD=-80+dx:dx:80-dx;
hold on; plot(0,0,'k+','MarkerSize',5);
for jj=[10] %Draw stimuli circles
    th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
    plot(xunit, yunit,'--k','LineWidth',1);
end
set(h,'Xdata',XD,'Ydata',YD );
xlim([-60+dx, 60-dx]);
ylim([-60+dx, 60-dx]);
c=colorbar;
% c.Label.String = ;
title(['initial \DeltaRF dist. from ',dRF_V1, ' for ', LM_Group2_RFs],'FontSize',8);
axis square


ax2=subplot(2,xPlots,2)
h=imagesc((counts_G1'));
set(gca,'YDir','normal')
dx=2.5; 
XD=-80+dx:dx:80-dx;
YD=-80+dx:dx:80-dx;
hold on; plot(0,0,'k+','MarkerSize',5);
for jj=[10] %Draw circles
    th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
    plot(xunit, yunit,'--k','LineWidth',1);
end
set(h,'Xdata',XD,'Ydata',YD );
xlim([-60+dx, 60-dx]);
ylim([-60+dx, 60-dx]);
c=colorbar;
title(['initial \DeltaRF dist. from ',dRF_LM, ' for ', LM_Group1_RFs],'FontSize',8);
axis square

% initial 2D histogram dRF for Group 2
ax8=subplot(2,xPlots,xPlots+2)
h=imagesc((counts_G2'));
set(gca,'YDir','normal')
dx=2.5; 
XD=-80+dx:dx:80-dx;
YD=-80+dx:dx:80-dx;
hold on; plot(0,0,'k+','MarkerSize',5);
for jj=[10] %Draw stimuli circles
    th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
    plot(xunit, yunit,'--k','LineWidth',1);
end
set(h,'Xdata',XD,'Ydata',YD );
xlim([-60+dx, 60-dx]);
ylim([-60+dx, 60-dx]);
c=colorbar;
% c.Label.String = ;
title(['initial \DeltaRF dist. from ',dRF_LM, ' for ', LM_Group2_RFs],'FontSize',8);
axis square

% 2D histogram survival dRF Group 2
ax4=subplot(2,xPlots,4)
h=imagesc((survivalG1));
set(gca,'YDir','normal')

dx=2.5; %160/64
XD=-80+dx:dx:80-dx;
YD=-80+dx:dx:80-dx;

set(h,'Xdata',XD,'Ydata',YD );
xlim([-60+dx, 60-dx]);
ylim([-60+dx, 60-dx]);
hold on; plot(0,0,'k+','MarkerSize',5);
c=colorbar;
for jj=[10] %Draw stimuli circles
    th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
    plot(xunit, yunit,'--k','LineWidth',1);
end
title(['Surviving \DeltaRF dist. - ' LM_Group1_RFs  ],'FontSize',8);
axis square

% 2D histogram survival dRF Group 2
ax10=subplot(2,xPlots,xPlots+4)
h=imagesc((survivalG2));
set(gca,'YDir','normal')

dx=2.5; %160/64
XD=-80+dx:dx:80-dx;
YD=-80+dx:dx:80-dx;

set(h,'Xdata',XD,'Ydata',YD );
xlim([-60+dx, 60-dx]);
ylim([-60+dx, 60-dx]);
c=colorbar;
hold on; plot(0,0,'k+','MarkerSize',5);
for jj=[10] %Draw stimuli circles
    th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
    plot(xunit, yunit,'--k','LineWidth',1);
end
hold off
title(['Surviving \DeltaRF dist. - ' LM_Group2_RFs  ],'FontSize',8);
axis square


% Fraction of RF overlap
ax3=subplot(2,xPlots,3)
binMeans_Group1_plot = binMeans_Group1;
binMeans_Group1_plot(isnan(binMeans_Group1_plot)) = 0;

h=imagesc(binMeans_Group1_plot);
set(gca,'YDir','normal')
dx=2.5; 
XD=-80+dx:dx:80-dx;
YD=-80+dx:dx:80-dx;
hold on; plot(0,0,'k+','MarkerSize',5);
for jj=[10] %Draw stimuli circles
    th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
    plot(xunit, yunit,'--k','LineWidth',1);
end
set(h,'Xdata',XD,'Ydata',YD );
xlim([-60+dx, 60-dx]);
ylim([-60+dx, 60-dx]);
c=colorbar;
c.Label.String = ['Fraction RF overlap'];
title({['Mean Fraction of RF Overlap in \DeltaRF space; '];[ LM_Group1_RFs]},'FontSize',8);
axis square

% Fraction of RF overlap
ax9=subplot(2,xPlots,xPlots+3)
binMeans_Group2_plot = binMeans_Group2;
binMeans_Group2_plot(isnan(binMeans_Group2_plot)) = 0;
h=imagesc(binMeans_Group2_plot);
set(gca,'YDir','normal')
dx=2.5;
XD=-80+dx:dx:80-dx;
YD=-80+dx:dx:80-dx;
hold on; plot(0,0,'k+','MarkerSize',5);
for jj=[10] %Draw stimuli circles
    th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
    plot(xunit, yunit,'--k','LineWidth',1);
end
set(h,'Xdata',XD,'Ydata',YD );
xlim([-60+dx, 60-dx]);
ylim([-60+dx, 60-dx]);
c=colorbar;
c.Label.String = ['Fraction RF overlap'];
title({['Mean Fraction of RF Overlap in \DeltaRF space; '];[LM_Group2_RFs]},'FontSize',8);
axis square

% dRF 2D histogram diference (G1-G2) fraction of overlap of RF per dRF location
ax5=subplot(2,xPlots,5)
h=imagesc(binMeans_Group1_plot-binMeans_Group2_plot,[-.1 .1]);
set(gca,'YDir','normal')
dx=2.5;
XD=-80+dx:dx:80-dx;
YD=-80+dx:dx:80-dx;
hold on; plot(0,0,'k+','MarkerSize',5);
for jj=[10] %Draw stimuli circles
    th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
    plot(xunit, yunit,'--k','LineWidth',1);
end
set(h,'Xdata',XD,'Ydata',YD );
xlim([-60+dx, 60-dx]);
ylim([-60+dx, 60-dx]);
colormap(divColorMap)
c=colorbar;
c.Label.String = ['\Delta Overlap (Vert-Horz)'];
title({['Mean Fraction of RF Overlap difference in \DeltaRF space'] ;[LM_Group1_RFs ' - ' LM_Group2_RFs  ]},'FontSize',8);
axis square

% dRF 2D histogram diference between Horizontal and Vertical Preferring
% surviving populations
ax11=subplot(2,xPlots,xPlots+5)
currMaxValue=max(abs(survivalG1/max(survivalG1,[],'all')-survivalG2/max(survivalG2,[],'all')),[],'all');
h=imagesc(survivalG1/max(survivalG1,[],'all')-survivalG2/max(survivalG2,[],'all'),[-currMaxValue currMaxValue]);
set(gca,'YDir','normal')

dx=2.5; 
XD=-80+dx:dx:80-dx;
YD=-80+dx:dx:80-dx;
hold on; plot(0,0,'k+','MarkerSize',5);
for jj=[10] %Draw stimuli circles
    th = 0:pi/50:2*pi; xunit = jj * cos(th) + 0; yunit = jj * sin(th) + 0;
    plot(xunit, yunit,'--k','LineWidth',1);
end
set(h,'Xdata',XD,'Ydata',YD );
xlim([-60+dx, 60-dx]);
ylim([-60+dx, 60-dx]);
c=colorbar;
c.Label.String = '\Delta Fraction of inputs (Vert-Horz)';
title({ ['Surviving \DeltaRF distribution'];[LM_Group1_RFs ' - ' LM_Group2_RFs  ]},'FontSize',8);
axis square


% dRF amp histogram
subplot(2,xPlots,xPlots)
clear h0 h1 dRF_amp_counts_init_G12 dRF_amp_counts_surv_G12 dRF_amp_counts_bins

[~, rhoG12]=cart2pol([rotated_samples_LM_Group1(:,1);rotated_samples_LM_Group2(:,1)],[rotated_samples_LM_Group1(:,2);rotated_samples_LM_Group2(:,2)]);
[~, rhoSurvivingG12]=cart2pol([rotated_samples_LM_Group1(logical(survival_group1'),1);rotated_samples_LM_Group2(logical(survival_group2'),1)] ,[rotated_samples_LM_Group1(logical(survival_group1'),2); rotated_samples_LM_Group2(logical(survival_group2'),2)]);
[dRF_amp_counts_init_G12(1,:),dRF_amp_counts_bins] = histcounts( rhoG12, 0:1:60,'normalization','probability');
[dRF_amp_counts_surv_G12(1,:),dRF_amp_counts_bins] = histcounts( rhoSurvivingG12 , 0:1:60,'normalization','probability');

hold on
h0=plot(dRF_amp_counts_bins(1:end-1),dRF_amp_counts_init_G12(1,:),'-','Color',init_color,'LineWidth',3);
h1=plot(dRF_amp_counts_bins(1:end-1),dRF_amp_counts_surv_G12(1,:),'-','Color','k','LineWidth',3);
xlim([dRF_amp_counts_bins(1) dRF_amp_counts_bins(end)])
ylabel('Fraction')
xlabel('||\DeltaRF||(º)')
set(gca,'xtick',0:10:100,'xlim',[0 60],'ylim',[0 0.07],'ytick',0:0.01:1);
legend([h0,h1],{'init','surv'})
title('Average of G1 and G2','FontSize',8)
hold off

% collinear depletion index
subplot(2,xPlots,xPlots+xPlots)
clear h1 h2
[h1,vert_data_surv]=plotOverlapTwoGroups_RDLP(survivalG2,survivalG1,10,0,'k',LM_Group2_RFs,LM_Group1_RFs);
hold on
[h2,vert_data_init]=plotOverlapTwoGroups_RDLP(counts_G2',counts_G1',10,0,init_color,LM_Group2_RFs,LM_Group1_RFs);
plot([0 5],[.5 .5],'k--')
set(gca,'ylim',[0 .75],'ytick',0:0.10:1)
legend([h1,h2],{['after survival'],['initial dist']},'location','NorthWest')
hold off

colormap(ax1,mainColorMap)
colormap(ax2,mainColorMap)
colormap(ax3,mainColorMap)
colormap(ax4,mainColorMap)
colormap(ax7,mainColorMap)
colormap(ax8,mainColorMap)
colormap(ax9,mainColorMap)
colormap(ax10,mainColorMap)

sgtitle({titleText; ['\DeltaRF starting distribution: ',dRF_LM,'; ',LM_Group1_RFs,' vs ',LM_Group2_RFs,'; ', connectivityRule, ' connectivity rule to V1; survival rate exponential factor: ',num2str(survivalGain), '; fraction of LMV1 joint RF']})
set(gcf,'Position', [10         10        1940         920]);

saveas(gcf,[saveFigFolder filesep 'simul_init_',dRF_LM,'_',dRF_V1,'_G1_',LM_Group1_RFs,'_G2_',LM_Group2_RFs,'_',V1_RFs,'_dRF_', num2str(numberOfSamples/1000),'k_',customShortText,'_ecc_',num2str(ecc_factor*10)],'fig')
export_fig([saveFigFolder filesep 'simul_init_',dRF_LM,'_',dRF_V1,'_G1_',LM_Group1_RFs,'_G2_',LM_Group2_RFs,'_',V1_RFs,'_dRF_', num2str(numberOfSamples/1000),'k_',customShortText,'_ecc_',num2str(ecc_factor*10)], '-png')


[~, rho_G1]=cart2pol(rotated_samples_LM_Group1(:,1),rotated_samples_LM_Group1(:,2));
[~, rho_G2]=cart2pol(rotated_samples_LM_Group2(:,1),rotated_samples_LM_Group2(:,2));

med_ecc = nanmedian([ecc_G1(logical(survival_group1)& rho_G1'>=10) ; ecc_G2(logical(survival_group2)& rho_G2'>=10)]);
mean_ecc = nanmean([ecc_G1(logical(survival_group1)& rho_G1'>=10) ; ecc_G2(logical(survival_group2)& rho_G2'>=10)]);

%% Save data 
save([saveFigFolder filesep 'simul_init_',dRF_LM,'_',dRF_V1,'_G1_',LM_Group1_RFs,'_G2_',LM_Group2_RFs,'_',V1_RFs,'_', num2str(numberOfSamples/1000),'k_',customShortText,'_ecc_',num2str(ecc_factor*10),'_data.mat'],...
    'vert_data_surv','vert_data_init','med_ecc','mean_ecc' )

end

%%
%% Supplementary Functions

%% Get parameters to generate 2D dRF distributions
function [AziSD, EleSD, twoDimGaussianAngle]=getdRFGaussianParameters(ExperimentName)
switch ExperimentName
    case 'DRP0'
        AziSD=10;
        EleSD=8;
        twoDimGaussianAngle=10;
    case 'V1'
        AziSD=6.7;
        EleSD=4.5;
        twoDimGaussianAngle=0;
end
end

%% Get parameters from experimental data to generate RF shape
function [RF_longAxis, RF_longAxisSD, RF_shortAxis, RF_shortAxisSD, RF_angle, RF_angleSD]=getRFParameters(RF_Group);

loadFolder = ['RFshape_parameters'];

if ~contains(RF_Group,'V1')
    if contains(RF_Group(1:6),'NR')
        if contains(RF_Group(1:6),'L5')
            load([loadFolder filesep 'RFshape_NR_L5.mat'])
        elseif contains(RF_Group(1:6),'L23')
            load([loadFolder filesep 'RFshape_NR_L23.mat'])
        else
            load([loadFolder filesep 'RFshape_NR.mat'])
        end

    elseif contains(RF_Group(1:6),'DRP21')
        if contains(RF_Group(1:6),'L5')
            load([loadFolder filesep 'RFshape_DRP21_L5.mat'])
        elseif contains(RF_Group(1:6),'L23')
            load([loadFolder filesep 'RFshape_DRP21_L23.mat'])
        else
            load([loadFolder filesep 'RFshape_DRP21.mat'])
        end

    elseif contains(RF_Group(1:6),'DRP0')
        load([loadFolder filesep 'RFshape_DRP0.mat'])
    end
end

if contains(RF_Group,'OS90')

    RF_longAxis    = RFshape.longAxis.mean.Or90;
    RF_longAxisSD  = RFshape.longAxis.std.Or90;
    RF_shortAxis   = RFshape.shortAxis.mean.Or90;
    RF_shortAxisSD = RFshape.shortAxis.std.Or90;
    RF_angle       = RFshape.rfangle.mean.Or90;
    RF_angleSD     = RFshape.rfangle.std.Or90;

elseif contains(RF_Group,'OS0')

    RF_longAxis    = RFshape.longAxis.mean.Or0;
    RF_longAxisSD  = RFshape.longAxis.std.Or0;
    RF_shortAxis   = RFshape.shortAxis.mean.Or0;
    RF_shortAxisSD = RFshape.shortAxis.std.Or0;
    RF_angle       = RFshape.rfangle.mean.Or0;
    RF_angleSD     = RFshape.rfangle.std.Or0;

elseif contains(RF_Group,'V1')
    RF_longAxis=16.1;
    RF_longAxisSD= 5.8;
    RF_shortAxis=6.9;
    RF_shortAxisSD=1.8;
    RF_angle=14.8;
    RF_angleSD=24.6;
end

end

%%
function overlap = calculateEllipseOverlap_LPRD(a1, b1, deg1, x1, y1, a2, b2, deg2, x2, y2, distance,borderOfEllipseFromSD,HalfVisualSpan)
%borderOfEllipseFromSD=1 means border considered for overlap is the SD of the RF gaussian fit.
TotalVisualSpan=HalfVisualSpan*2;
 
a1=a1*borderOfEllipseFromSD;
b1=b1*borderOfEllipseFromSD;

a2=a2*borderOfEllipseFromSD;
b2=b2*borderOfEllipseFromSD;

theta1=deg2rad(deg1);
theta2=deg2rad(deg2);

% Generate points on ellipses
t = linspace(0, 2*pi, 100);

ellipse1_x = x1 + a1 * cos(t) * cos(theta1) - b1 * sin(t) * sin(theta1);
ellipse1_y = y1 + a1 * cos(t) * sin(theta1) + b1 * sin(t) * cos(theta1);

ellipse2_x = x2 + a2 * cos(t) * cos(theta2) - b2 * sin(t) * sin(theta2);
ellipse2_y = y2 + a2 * cos(t) * sin(theta2) + b2 * sin(t) * cos(theta2);

% Translate one ellipse to the specified distance
ellipse2_x = ellipse2_x + distance;

E1=poly2mask(ellipse1_x,ellipse1_y,TotalVisualSpan,TotalVisualSpan);
E2=poly2mask(ellipse2_x,ellipse2_y,TotalVisualSpan,TotalVisualSpan);

SumE=E1+E2;
pixOverlap=sum(SumE==2,'all');
jointArea=SumE>0;

%fraction of ellipse overlap  over the joint RFs
overlap=pixOverlap/sum(jointArea,'all');

end

%%
function rotated_samples=sampleFromRotated2DGaussian_LP(mu,sigma_x,sigma_y,rotation_angle_in_deg,num_samples)
%% Parameters of the Gaussian distribution
% mu = [0, 0];               % Mean
% sigma_x = 2;               % Standard deviation along x-axis
% sigma_y = 1;               % Standard deviation along y-axis
% rotation_angle = 10;       % Rotation angle in degrees

% Converting angle to radians
rotation_angle=deg2rad(rotation_angle_in_deg);
% Covariance matrix
sigma = [sigma_x^2, 0; 0, sigma_y^2];
% Rotation matrix
R = [cos(rotation_angle), -sin(rotation_angle); sin(rotation_angle), cos(rotation_angle)];
% Sample from the 2D Gaussian distribution
samples = mvnrnd(mu, sigma, num_samples);
% Rotate the samples
rotated_samples = (R * samples')';

end

%%
function [ numberHorizontalBin, numberVerticalBin, numMatched]=getAngularAverages_dRF_Data_LP(binMeans,distantThresh, binAngle)

mat=binMeans;
cent = size(mat)./2;
[rr,cc] = meshgrid(1:size(mat,1),1:size(mat,2));
[TH,R,Z] = cart2pol(rr(:)-cent(1),cc(:)-cent(2),mat(:));
polarCoordinates = [TH R Z];

matTHdeg=reshape(rad2deg(TH),size(rr));
matR=reshape(R,size(rr));

verticalPositions=(matTHdeg>-135+binAngle&matTHdeg<-45+binAngle)|(matTHdeg>45+binAngle&matTHdeg<135+binAngle);
verticalAndDistantPositions=verticalPositions&matR>=distantThresh;
horizontalAndDistantPositions=~verticalPositions&matR>=distantThresh;

numberVerticalBin=sum(binMeans(verticalAndDistantPositions));
numberHorizontalBin=sum(binMeans(horizontalAndDistantPositions));
numMatched=sum(binMeans(matR<distantThresh));

end

%%
function [h, dataPoints]=plotOverlapTwoGroups_RDLP(binMeans1, binMeans2, distantThresh,binAngle,col,g1,g2)

[H1, V1, M1]=getAngularAverages_dRF_Data_LP(binMeans1, distantThresh,binAngle);
[H2, V2, M2]=getAngularAverages_dRF_Data_LP(binMeans2, distantThresh,binAngle);


h=plot([1 2],[V1/(H1+V1) , V2/(H2+V2)],'color',col)
ax=gca;
ylim(ax,[0.15 0.6])
xlim(ax,[0.9 2.1])
set(ax,'Xtick',1:2,'XtickLabel',{g1, g2})
ylabel(ax, 'Fraction of boutons deviating in elevation')

dataPoints = [V1/(H1+V1) , V2/(H2+V2)];

end