clear all

%% TEST PARAMETERS

%TESTING NANOVNA Code, BSA and antibsa test, 0.01% anti-BSA
file_location = "CHANGE ME";
num_samps = 1024;
determine_binding_coeffs = 1;
drift_correction = 1;
molar_conc = 1e-6; %concenration of protein in M/L
index_dissoc_start = 500; index_dissoc_end = 1000; index_assoc_start = 10; index_assoc_end = 400; %These are the data indicies associated with the start of the association and dissociation phases of the protein binding
drift_correct_index_start = 1; drift_correct_index_end = 10; %these are indicies used to correct for any drift in the response. These should be chosen from a portion of the frequency over time that is somewhat flat

%% Analysis
files = dir(file_location);

%initializing variables
data1 = zeros(num_samps , 3 , length(files)-2);
s11 = zeros(num_samps , length(files)-2);
Z1 = zeros(length(files)-2 , num_samps);
q_fac_1 = zeros(num_samps , length(files)-2);
day = zeros(1 , length(files)-2);
hour = zeros(1 , length(files)-2);
minute = zeros(1 , length(files)-2);
second = zeros(1 , length(files)-2);
min_val1 = zeros(1,length(files)-2);
min_index1 = zeros(1,length(files)-2);
zero_freq1 = zeros(1,length(files)-2);
centroid_1 = zeros(1 , length(files)-2);
min_freq1 = zeros(1,length(files)-2);
first_res_mag_1 = zeros(1,length(files)-2);

for index = 1:(length(files)-2)
    file_names{index} = files(index+2).name;
    [data1(:,:,index) , s11(:,index) , Z1(index,:) , q_fac_1(:,index) , centroid_1(:,index)] = process_data(strcat(file_location,file_names{index}));

    [~, ~, day(:,index), hour(:,index), minute(:,index), second(:,index)] = s1p_name_to_datetime(file_names{index});

    [min_val1(index) , min_index1(index)] = min(s11(:,index));

    %simple and fitting_algorithm are the two possible arguments. Simple takes less time, but may include more noise wheras fitting_algorithm uses a linear regression algorithm to more accurately determine the zero point. simple is ~ 30 us per data point, fitting_algorithm is ~ 10 ms
    [zero_freq1(index)] = find_zero_phase(data1(:,1,index) , data1(:,3,index) , "simple"); 
    
    min_freq1(index) = data1(min_index1(index),1,index);

    if index ~= 1
        first_res_mag_1(index) = s11(first_min_freq_1,index);
    else
        [first_res_mag_1(1) , first_min_freq_1] = min(s11(:,1));
    end
end

diff_second_isnan = isnan(diff(second));

diff_second_isnan_location = find(diff_second_isnan);

diff_second = diff(second);

diff_second(diff_second_isnan_location) = 0;

increment = median(diff_second)/60;

time = (0:index-1)*increment;

figure
plot((0:index-1)*increment , (centroid_1 - centroid_1(1,1)));
grid on
hold on
xlabel('Time (min)');
ylabel('Resonance frequency (MHz)');
title('Centroid Shift of Reflection Coefficent Over Time');
set(gcf, 'Units','Inches','Position',[1,1,9,5.75],'PaperUnits','Inches', 'PaperSize',[9,5.75]);

figure
plot((0:index-1)*increment,zero_freq1)
hold on
grid on
xlabel('Time (min)')
ylabel('Resonance Frequency (MHz)')
title('Resonance Shift determined from Phase Over Time');
set(gcf, 'Units','Inches','Position',[1,1,9,5.75],'PaperUnits','Inches', 'PaperSize',[9,5.75]);

figure
plot(zero_freq1)
grid on
xlabel('Data Index');
ylabel('Resonance Frequency (MHz)')
title('Resonance Shift determined from Phase with Data Index');
set(gcf, 'Units','Inches','Position',[1,1,9,5.75],'PaperUnits','Inches', 'PaperSize',[9,5.75]);


%This determines the fits for the binding kinetics based on chosen indicies
%from the collected data
if determine_binding_coeffs == 1
    if drift_correction == 1
        binding_curve = zero_freq1 - (zero_freq1(drift_correct_index_start) - zero_freq1(drift_correct_index_end)) / (time(drift_correct_index_start) - time(drift_correct_index_end)) * time;
    else
        binding_curve = zero_freq1;
    end
    [fittest_diss,gof_diss] = calc_dissociation(index_dissoc_start, index_dissoc_end , index_assoc_start,60* time, binding_curve, 1);
    [fittest_assoc,gof_assoc] = calc_association(index_assoc_start, index_assoc_end, 60*time, binding_curve, molar_conc, fittest_diss.b, 1);
    plot_kinetic_data(time * 60 , binding_curve,...
        fittest_assoc , fittest_diss , index_assoc_start , index_assoc_end,...
        index_dissoc_start , index_dissoc_end)
    xlabel('Time (s)');
    ylabel('Relative Resonance Shift');
    title('Resonance Shift with Overlaid Binding Fits')
    Kd = fittest_diss.b / (fittest_assoc.b*1e5);
    disp(['Kd = ', num2str(Kd * 1e9),' nM']);
end


function [centroid_freq] = centroid_func(centroid_indicies, data, frequency)
centroid_freq = sum(data(centroid_indicies) .* frequency(centroid_indicies))./sum(data(centroid_indicies));
end

function [data , s11_test , Z , q_fac , centroid] = process_data(file_location)

data = readmatrix(file_location , 'FileType','text');

data(:,2) = 20*log10(data(:,2));

%Generating Impedance from S11
s11_test = data(:,2);


[real_refl,imag_refl] = pol2cart(data(:,3)*pi/180 , 10.^(data(:,2)./20));

complex_refl = (real_refl + 1i*imag_refl);

Z = 50*(1 + complex_refl)./(1 - complex_refl);

%FWHM
[~ , min_ind ] = min(s11_test);
centroid_indicies = find((data(:,1) >= data(min_ind , 1) - 10e6) & (data(:,1) <= data(min_ind , 1) + 10e6)); %Sums values withint 20 MHz of peak

%need to sum inverted values, this assumes they are reflectivity data
centroid = centroid_func(centroid_indicies , 1-10.^(data(:,2)./20) , data(:,1));

[peak_val , peak_index] = min(10.^(data(:,2)./20));

% FWHM_indicies = find((10.^(data(:,2)./20) <= 0.5*(max(10.^(data(:,2)./20)) - peak_val)));
FWHM_indicies = find((10.^(data(:,2)./20) <= 0.5*(2*peak_val - 10.^(data(:,2)./20) + max(10.^(data(:,2)./20)))));

q_fac = data(peak_index,1)./(data(FWHM_indicies(end),1) - data(FWHM_indicies(1),1));
end

function [S11_peak_mag , S11_peak_freq , admit_peak_mag , admit_peak_freq] = calc_peak_freq(data , in_dB , in_impedance)

if in_impedance == 1
    Z = data(:,2);
    s11_test = (data(:,2).*(cosd(data(:,3)) + 1i*sind(data(:,3))) - 50)./(data(:,2).*(cosd(data(:,3)) + 1i*sind(data(:,3))) + 50);
else
    s11_test = data(:,2);
    s11_ang = data(:,3);

    if in_dB == 1
        [real_refl,imag_refl] = pol2cart(s11_ang*pi/180 , 10.^(s11_test./20));
    else
        [real_refl,imag_refl] = pol2cart(s11_ang*pi/180 , s11_test);
    end

    complex_refl = real_refl + 1i*imag_refl;

    Z = 50*(1 + complex_refl)./(1 - complex_refl);
end



[S11_peak_mag , s11_peak_ind] = min(s11_test);
[admit_peak_mag, admit_peak_ind] = max(abs(1./Z));

S11_peak_freq = data(s11_peak_ind,1);
admit_peak_freq = data(admit_peak_ind,1);
end