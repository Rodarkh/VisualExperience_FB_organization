function generate_Figure_6A(initial_LM_dRF, initial_V1_dRF, dendriticSize, V1_RFshape, connectivityRule, borderOfEllipseFromSD, survivalGain, nSamples,ec, nIter,motherSaveFolder) 
% Model created and implemented by Rodrigo F. Dias and Leopoldo Petreanu
% Code written by Rodrigo F. Dias and Leopoldo Petreanu - Cortical Circuits Laboratory - Champalimaud Research - 2024
%
% Work under revision (NEURON-D-22-01431R1)
% "Visual experience reduces the spatial redundancy between cortical feedback inputs and primary visual cortex neurons." 
%
%% settings
currSettings = ['dend ' ,num2str(dendriticSize), '; surv ',num2str(survivalGain), '; ',num2str(connectivityRule), '; SD ',num2str(borderOfEllipseFromSD)];

groups={'NR_L5','NR'};
color2plot={rgb('Navy');[0 0 0]};
nGroups = length(groups);

clear listing_save data2load ecc_median ecc_mean 
dataFolder=[motherSaveFolder, filesep, 'Model', filesep, currSettings, '_AspectRatio'];  
cd(dataFolder);

listing_save = dir(dataFolder);

for iii=1:nGroups
    for iter=1:nIter
        for j=1:length(ec)
            clear currFile loadedData
            currFile = ['simul_init_',initial_LM_dRF,'_',initial_V1_dRF,'_G1_',groups{iii},'_OS90_G2_',groups{iii},'_OS0_',V1_RFshape,'_',num2str(nSamples/1000),'k__iter',num2str(iter),'_ecc_',num2str(ec(j)*10),'_data.mat'];

            loadedData = load([dataFolder filesep currFile]);

            ecc_median(iii,j,iter) = loadedData.med_ecc;
            ecc_mean(iii,j,iter) = loadedData.mean_ecc;
            ver_data_surv(iii,:,j,iter) = loadedData.vert_data_surv;
            ver_data_init(iii,:,j,iter) = loadedData.vert_data_init;
        end
    end
end
ver_dif_surv=squeeze([ver_data_surv(:,1,:,:) - ver_data_surv(:,2,:,:)]);
ver_dif_init=squeeze([ver_data_init(:,1,:,:) - ver_data_init(:,2,:,:)]);

ar_median = sqrt( 1 ./ (1-ecc_median.^2) );

%%
cd ..

for i=1:nGroups
    for iter=1:nIter
        for j=1:length(ec)
            currData{i}(iter,j) = ver_dif_surv(i,j,iter);
        end
    end
end

figure('Position', [0 40 900 920]);
hold on;
for i=1:nGroups
    plotLWEA(nanmean(ar_median(i,:,:),3),nanmean(ver_dif_surv(i,:,:),3),nanstd(ver_dif_surv(i,:,:),0,3),color2plot{i},.3,[])
end
for i=1:nGroups
    plot(nanmean(ar_median(i,:,:),3),nanmean(ver_dif_surv(i,:,:),3),'-','color',color2plot{i},'LineWidth',3)
    plot(nanmean(ar_median(i,:,:),3),nanmean(ver_dif_surv(i,:,:),3),'o','color',color2plot{i},'LineWidth',2,'MarkerSize',8)
    plot(nanmean(ar_median(i,3,:),3),nanmean(ver_dif_surv(i,3,:),3),'o','color',color2plot{i},'LineWidth',2,'MarkerSize',8,'markerfacecolor',color2plot{i})
end
plot([0 6],[0 0],'k--')
ylabel(['Collinear depletion index'])
xlabel('Median aspect ratio')
title('Model')
set(findobj(gcf,'type','axes'),'FontSize',16,'xlim',[.99 4.1],'XTick',0:.5:5,'ylim',[-.015 .081],'ytick',-0.2:0.01:2);

saveas(gcf,['Figure_6A'], 'fig')
saveas(gcf,['Figure_6A'], 'svg')
export_fig(['Figure_6A'], '-png')


end


