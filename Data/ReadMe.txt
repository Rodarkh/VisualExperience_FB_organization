Visual experience reduces the spatial redundancy between cortical feedback inputs and primary visual cortex neurons.
Rodrigo F. Dias*, Radhika Rajan*, Margarida Baeta, Beatriz Belbut, Tiago Marques and Leopoldo Petreanu†.
Champalimaud Neuroscience Programme, Champalimaud Foundation, Lisbon, Portugal. Neuron, 2024


This repository contains the main dataset used in this paper. 

Each data structure contains data for more than one group:

DRP0_DRP21_NR_dataset.mat - 3 wild type animal groups with different visual experiences;
L23_DRP21_NR_dataset.mat  - 2 cux2 (feedback originating from LM L2/3) animal groups with different visual experiences;
L5_DRP21_NR_dataset.mat   - 2 Rbp4 (feedback originating from LM L5) animal groups with different visual experiences;

Most variables are organized as cells, with the first dimension of the cell representing the visual experience, the different groups. 
Each structure has variables with the same names. 
Unless stated in the variable name (by including "V1"/"N"/"neur"), each will have data per BOUTON or per ANIMAL. 
Here we detail the contents of these data structures:

currGroups = name of the groups of current structure. Identifies order of groups inside variables.
boot_RF_CI_xy = Receptive Field 95% Confidence Interval distance in x and y coordinates
boot_RF_CI_amp = Receptive Field 95% Confidence Interval amplitude (euclidean distance)
meanN = mean V1 neurons Receptive Field center per imaging session
neuronFitData = matrix with data per neuron. Includes the 9 parameters resulting from the 2D gaussian fit to estimate the receptive field. (offset , amplitude, xCenter, yCenter, RotationAngle, xAxis, yAxis, sum of squares, Rsquare)
nSess = ID number of session, resets for each animal
nSession = total number of session per animal
nTotalROIs_ani = total number of boutons per animal before calculating Receptive Field
real_animalID = animal ID number per bouton
real_animalID_neur_all = animal ID number per neuron (all neuron ROIs, regardless of Receptive Field quality)
real_boutonFitResults_all = matrix with data for all bouton ROIs. Includes the 9 parameters resulting from the 2D gaussian fit to estimate the receptive field. (offset , amplitude, xCenter, yCenter, RotationAngle, xAxis, yAxis, sum of squares, Rsquare)
real_dRF = individual Receptive Field center per bouton minus mean Receptive Field center for neurons of same session (aximuth and elevation coordinates)
real_dRF_amp = individual Receptive Field center per bouton minus mean Receptive Field center for neurons of same session (euclidian/radial distance)
real_dRF_tpos = individual Receptive Field center per bouton minus mean Receptive Field center for neurons of same session (angular difference where 0 is 0 azimuth and along the elevation axis, increasing clockwise)
real_DSpind= logical vector, true if Direction Selective (t-test significant at 5% of preferred vs anti-preferred)
real_DSI = Direction Selectivity Index - calculated from biggest response
real_DSIPolar = Direction Selectivity Index - calculated from vectorial sum of the responses to the 8 directions. Normalized amplitude of the resulting angle is the DSIPolar
real_ec = eccentricity value (sqrt(1- short^2/longAxis^2) ) 
real_goodRF_all = logical values corresponding to Rsquare>0.35;
real_grnOSDS = logical values signaling grating responsive boutons that are not OS nor DS
real_grRespind = logical values signaling grating responsive boutons (SNR>1)
real_longRFwidth = value in degrees of the long Receptive Field axis
real_meanResp = mean dF/F value for each grating direction
real_MeanTraceDir = mean dF/F trace for each grating direction (frame rate 6Hz, stimulus starts at frame7)
real_neuronFitData_all = matrix with data per neuron for all ROIs. Includes the 9 parameters resulting from the 2D gaussian fit to estimate the receptive field. (offset , amplitude, xCenter, yCenter, RotationAngle, xAxis, yAxis, sum of squares, Rsquare)
real_OSIpolar = Orientation Selectivity Index: calculated from vectorial sum of the responses to the 4 orientation. Normalized amplitude of the resulting angle is the OSIPolar
real_OSpind = logical vector, true if Orientation Selective (t-test significant at 5% of preferred vs orthogonal)
real_prefDir = index of preferred direction (1 is upwards, grows clockwise)
real_prefDirPolar = preferred direction calculated from vectorial sum of the responses to the 8 directions.
real_prefOrPolar = preferred orientation calculated from vectorial sum of the responses to the 4 orientations.
real_prefOrPolar_neur_all = preferred orientation calculated from vectorial sum of the responses to the 4 orientations.
real_RFangle_corr = Receptive Field angle with positive azimuth angle (grows clockwise)
real_RFarea = Receptive Field area
real_RS = Rsquare from Receptive Field fit
real_sessionID = Session ID per bouton, resets per animal
real_shortRFwidth = value in degrees of the short Receptive Field axis
real_SNR = SNR value for grating responses
real_x0_B = azimuth coordinate (degrees) for RF center
real_y0_B = elevation coordinate (degrees) for RF center
totalROIs_all = total number ROIs per group
