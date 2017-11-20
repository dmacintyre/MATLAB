clear all;
clc;
close all;

for phase = -1:0.05:1

clf;
N = 32;
SAMPLING_PHASE = phase;
RESAMPLE_FILTER_TAPS = 2000;
RESAMPLE_FILTER_BETA = 1;

expected = [50,50,50,50,50,50,50,47,35,16,6,16,35,47,50,50,50,50,50,50,50,47,35,16,3,1,1,1,1,1,-3,-16,-35,-47,-50,-47,-35,-16,-6,-15,-37,-50,-35,1,33,47,50,50,50,47,33,1 ,-33,-47,-50,-50,-50,-47,-33,1,35,50,37,16];
waveform = repmat(expected,1,30);
x = waveform;
n = 1:length(x);

% Generate ideal sample stream
% Running it through the resampler to maintain the same group delay through
% the filter
ideal_stream = non_int_resample(x,0,RESAMPLE_FILTER_BETA,RESAMPLE_FILTER_TAPS);

% Simulate a phase change due to a non-ideal sampling phase
resampled_stream = non_int_resample(x,SAMPLING_PHASE,RESAMPLE_FILTER_BETA,RESAMPLE_FILTER_TAPS);

% stem(ideal_stream);
% hold on
% stem(resampled_stream)

% At what rate are we going around the circle (radians per sample)
phase_change = (2*pi) / N;

% Wrap the input signal around the unit circle at a rate of N To determine
% magnitude and amplitude of N rate cos that is present
current_phase = phase_change*n;
val = ideal_stream.*exp(1j*current_phase);
ideal_calculated_phase = angle(sum(val));

% Wrap the resampled input  signal around the unit circle at a rate of N To determine
% magnitude and amplitude of N rate cos that is present
res_val = resampled_stream.*exp(1j*current_phase);
resampled_calculated_phase = angle(sum(res_val));

% Determine phase offset between ideal stream and resampled signal
resample_phase_error = ideal_calculated_phase - resampled_calculated_phase;

% Convert the phase offset into a number of samples delay
% This should match the value of SAMPLING_PHASE
samples_error = resample_phase_error / phase_change;

% To reconstruct the original waveform again filter with non-integer delay
% need to account for group delay again.

y_recover = non_int_resample(resampled_stream,samples_error,RESAMPLE_FILTER_BETA,RESAMPLE_FILTER_TAPS);
y_ideal = non_int_resample(ideal_stream,0,RESAMPLE_FILTER_BETA,RESAMPLE_FILTER_TAPS);

% Calculate Error
recovered_error = sum(abs(y_recover-y_ideal))./length(n)
unrecovered_error = sum(abs(resampled_stream-ideal_stream))./length(n)
improvement = 20*log10(unrecovered_error/recovered_error)

subplot(2,1,1)
stem(ideal_stream(100:200))
hold on
stem(resampled_stream(100:200))
legend('Ideal', 'Sampled')
title('No Timing Error Correction')

subplot(2,1,2)
stem(y_ideal(100:200))
hold on
stem(y_recover(100:200))
legend('Ideal', 'Recovered')
title('With Timing Error Correction')

error = sum(abs(y_recover-y_ideal))./length(n);
pause(1);

end;
