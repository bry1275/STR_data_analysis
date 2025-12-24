# STR_data_analysis
This code is used to extract data from a nanoVNA and plot the resonance frequency of an STR sensor. In addition, it also has the ability to correct for any drift over time as well as determine binding kinetics of any potential protein under study. This was developed by Bryan Chantigian at the University of Minnesota - Twin Cities for use in protein binding kinetic analysis in conjunction with an STR sensor as an electronic SPR analog. 

# STR_analysis_script
This is the main file used for data analysis and resonance frequency tracking. In order to use this code, the file location for where the nanoVNA code is stored must be put into the "file_location" variable. The number of VNA samples must also be added in the "num_samps" variable. If this is just used to track the progress of the experiment, the variables "determine_binding_coeffs" and "drift_correction" can be set to 0, as these are just flags for determining binding kinetics. To track the resonance shift over time, this script needs to be run.

For determining binding kinetics, the indicies in which the association phase start and stops, as well as the dissociation phase start and stop must be noted. These can be found from the plot "Resonance Shift determined from Phase with Data Index".

If any drift is noticed and impacts the binding kinetics, the variables "drift_correction_index_start" and "drift_correction_index_end" can be added. These should be found at a point in the test where the solution isn't changed and the drift is relatively linear. This will then correct for the drift when determining the binding kinetics.

Future iterations will constantly read the data and determine resonance frequency.
