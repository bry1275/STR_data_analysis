function [fittest_diss,gof_diss] = calc_dissociation(start_index, end_index, start_associatin_index, time, diss_curve, number_of_averages)
%CALC_DISSOCIATION is a function that determines that fits a dissociation
%curve to data.

%defining model parameters and equation
a = movmean(diss_curve(start_index)' - diss_curve(end_index), number_of_averages)*2; %Defining the coefficient this way was found to lead to better fits
c = movmean(diss_curve(end_index)' - diss_curve(start_associatin_index), number_of_averages)-a/2; %Defining the coefficient this way was found to lead to better fits
fit_model = fittype([num2str(a),'*exp(-b*x)+',num2str(c)]);

%fitting the model
[fittest_diss ,...
    gof_diss] = fit( time(start_index:end_index)' - time(start_index),...
    movmean(diss_curve(start_index:end_index)' - diss_curve(start_associatin_index), number_of_averages) ,...
    fit_model,'Lower',[0 0]);


end

