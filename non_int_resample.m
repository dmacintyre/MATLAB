function [ y ] = non_int_resample( x, delta, beta, N )
% x - Input signal to resample
% delta - Non-Integer amount to resample by
% beta - Used for Kaiser windowing
% N - number of filter taps

n = -N/2:1:N/2;
window = kaiser(length(n),beta);

h_resample = sinc(n-delta);
h = h_resample .* window';
y = conv(x,h,'same');

end

