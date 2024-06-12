% Model created and implemented by Rodrigo F. Dias and Leopoldo Petreanu
% Code written by Rodrigo F. Dias and Leopoldo Petreanu - Cortical Circuits Laboratory - Champalimaud Research - 2024
%
% Work under revision (NEURON-D-22-01431R1)
% "Visual experience reduces the spatial redundancy between cortical feedback inputs and primary visual cortex neurons." 
%

%% Run these settings to generate figures 5 and 6

clear
%% Model variables
n_iterations = 10;

initial_LM_dRF = 'DRP0';
initial_V1_dRF = 'V1';
dendriticSize = 15;
LM_RF_shape = {'DRP0','DRP21','NR','DRP21_L23','NR_L23','DRP21_L5','NR_L5'};
V1_RFshape = 'V1';
connectionRule =  'isoAngle';
SD = 2.5;
survivalGain = 15;
nSamples = 5000000;

saveFolder = 'Model Code';

%%  Run simulation loop
for iteration=1:n_iterations
    iteration
    for k=1:length(LM_RF_shape)
        Model_RFoverlap(initial_LM_dRF, initial_V1_dRF, dendriticSize, [LM_RF_shape{k},'_OS90'], [LM_RF_shape{k},'_OS0'], V1_RFshape, connectionRule, SD, survivalGain, nSamples, iteration,saveFolder)
        close all
    end
end

%% Generate Figure 5 using data from iterations

generate_Figure_5(initial_LM_dRF, initial_V1_dRF, dendriticSize, V1_RFshape, connectionRule, SD, survivalGain, nSamples, n_iterations,saveFolder)


%% Run simulation loop for Aspect Ratio modified Model

%using same conditions as above, but only running for the two groups compared in Figure 6

LM_RF_shape = {'NR_L5','NR'};
eccFactors = [0.6 0.8 1 1.2];

%%  Run simulation loop
for iteration=1:n_iterations
    iteration
    for k=1:length(LM_RF_shape)
        for iii=1:length(eccFactors)
            Model_RFoverlap_AspectRatio_mod(initial_LM_dRF, initial_V1_dRF, dendriticSize, [LM_RF_shape{k},'_OS90'], [LM_RF_shape{k},'_OS0'], V1_RFshape, connectionRule, SD, survivalGain, nSamples, eccFactors(iii), iteration,saveFolder)
        end
        close all
    end
end

%% Generate Figure 6

generate_Figure_6A(initial_LM_dRF, initial_V1_dRF, dendriticSize, V1_RFshape, connectionRule, SD, survivalGain, nSamples,eccFactors, n_iterations,saveFolder)
