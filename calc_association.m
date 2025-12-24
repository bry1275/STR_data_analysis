function [fittest_assoc,gof_assoc] = calc_association(start_index, end_index, time, assoc_curve, molar_concentration, diss_const, number_of_averages)
%CALC_ASSOCIATION is a function that determines that fits an association
%curve to data.

%defining the model equation
fit_model = fittype(['(a).*(1 - exp(-(b*1e5.*',num2str(molar_concentration),' + ',num2str(diss_const),').*(x)))']);
%Fitting the model
[fittest_assoc ,...
    gof_assoc] = fit( time(start_index:end_index)' - time(start_index),...
    movmean(assoc_curve(start_index:end_index)'-assoc_curve(start_index), number_of_averages) ,...
    fit_model, 'Lower',[0 0]);


end

