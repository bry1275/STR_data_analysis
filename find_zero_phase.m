function [zero_freq] = find_zero_phase(freq , phase_info , method)
%This function determines where the phase crosses zero using linear
%extrapolation.

if method == "simple"
    phase_more_zero = find(phase_info > 0);

    %doing linear extrapolation to find precise 0 crossing
    if isempty(phase_more_zero) ~= 1 && phase_more_zero(1) ~= 1
        slope = (phase_info(phase_more_zero(1)) - phase_info(phase_more_zero(1)-1))/(freq(phase_more_zero(1))-freq(phase_more_zero(1)-1));

        freq_change = - phase_info(phase_more_zero(1)-1)/slope;

        zero_freq = freq(phase_more_zero(1)-1) + freq_change;
    else
        zero_freq = NaN;
    end

elseif method == "fitting_algorithm"
    unwrapped_phase = 180/pi*unwrap(pi/180*phase_info);

    unwrapped_phase_less_zero = find(unwrapped_phase > 0);

    if isempty(unwrapped_phase_less_zero) == 1
        unwrapped_phase = unwrapped_phase + 180;
        unwrapped_phase_less_zero = find(unwrapped_phase < 0);
    end

    if isempty(unwrapped_phase_less_zero) == 1
        start_less_than_zero = int16.empty;
    else
        start_less_than_zero = unwrapped_phase_less_zero(1);
    end

    if length(start_less_than_zero) > 0 && start_less_than_zero(1) ~= 1 && start_less_than_zero > 50 && (start_less_than_zero < (length(unwrapped_phase)-50))
        [fittest , gof] = fit(freq(start_less_than_zero-50:start_less_than_zero+49) , ...
            unwrapped_phase(start_less_than_zero-50:start_less_than_zero+49),'poly1');

        zero_freq = -fittest.p2./fittest.p1;
    elseif length(start_less_than_zero) > 0 && start_less_than_zero(1) ~= 1 && start_less_than_zero < 50
        [fittest , gof] = fit(freq(start_less_than_zero-start_less_than_zero+1:start_less_than_zero+49) , ...
            unwrapped_phase(start_less_than_zero-start_less_than_zero+1:start_less_than_zero+49),'poly1');

        zero_freq = -fittest.p2./fittest.p1;
    elseif length(start_less_than_zero) > 0 && start_less_than_zero(1) ~= 1 && start_less_than_zero > length(unwrapped_phase)-50
        [fittest , gof] = fit(freq(start_less_than_zero-50:length(unwrapped_phase)) , ...
            unwrapped_phase(start_less_than_zero-50:length(unwrapped_phase)),'poly1');

        zero_freq = -fittest.p2./fittest.p1;
    else
        zero_freq = NaN;
    end

end
end