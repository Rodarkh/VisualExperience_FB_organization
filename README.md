Visual experience reduces the spatial redundancy between cortical feedback inputs and primary visual cortex neurons.
Rodrigo F. Dias*, Radhika Rajan*, Margarida Baeta, Beatriz Belbut, Tiago Marques and Leopoldo Petreanu†.
Champalimaud Neuroscience Programme, Champalimaud Foundation, Lisbon, Portugal


In this repository we have the Matlab functions necessary to run the model developed for this work, presented in Figure 5. 
Download all functions to a "VisualExperience_FB_organization-main" directory, and all data structures from the zip file to a "RFshape_parameter" directory inside "VisualExperience_FB_organization-main".
Changing Matlab to directory where functions are will allow to run Master_Script.mat and generate Figures 5 and 6 (takes several hours).


Contents:
Master_Script - this script will run the model with the conditions used to generate Figures 5 and 6. 

Model_RFoverlap - function that runs a single iteration of the model for the conditions specified in its inputs. 
It will get the parameters of the specified groups from the folder RFshape_parameters.

generate_Figure_5 - generates pannels C - J of Figure 5, using model data created by running Model_RFoverlap.

Model_RFoverlap_AspectRatio_mod - function that runs a single iteration of the model for the conditions specified in its inputs. 
It allows to change the aspect ratio of the Receptive Fields to study its effect. 

generate_Figure_6A -  generates pannels A of Figure 6, using model data created by running Model_RFoverlap_AspectRatio_mod.
