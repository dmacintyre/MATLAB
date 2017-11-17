% Project 2 - Q1
% Donald MacIntyre - djm4912

clear all;
close all;
clc;

% Define Filter Paramaters
fs = 20e3;
fp = 2e3;
fstop = 2.8e3;
ts = 1/fs;
deltap = 0.1088;    % -1 dB
deltas = 0.0032;    % -50 dB
wp = 2*pi*(fp/fs);
ws = 2*pi*(fstop/fs);

% Bilinear warping of discrete frequency omega = 2*tan(w/2)
warp_wp = 2*tan(wp/2);
warp_ws = 2*tan(ws/2);

% Calculate Passband ripple and stopband attenuation
Rp = -20*log10(1-deltap);
As = -20*log10(deltas);

% Calculate value of N necessary to meet the filter spec using Butterworth
top_n = (log10(10.^(Rp/10)-1)-log10(10.^(As/10)-1));
bottom_n = 2*log10(warp_wp/warp_ws);
N_pure = top_n / bottom_n;
N = ceil(N_pure);   % round to nearest integer

% Calculate IIR Butterworth Continous poles by solving 1+(s/j*cutoff)^2N=0
right = (((1/deltas).^2) - 1).^(1/(2*N));
warp_wc = warp_ws / right;
cutoff = 2*atan(warp_wc/2);
k = 0:1:2*N-1;
z = [];
p = 1j*warp_wc*exp(1j*(pi/(2*N))*(1+(2*k)));
p = p(real(p) < 0);     % Take just the poles of Hc(s)
k = real(prod(-p));     % Determine gain

[zd, pd, kd] = bilinear(z',p',k,1); % Bilinear mapping back to DT domain
[b,a] = zp2tf(zd,pd,kd);

[h,w] = freqz(b,a);

norm_w = w ./pi;
figure;plot(norm_w,abs(h));
grid on;
title('Filter Frequency Response')
xlabel('Normalized Frequency')
ylabel('Amplitude')
% figure;plot(norm_w,angle(h));
% figure;grpdelay(b,a)

% Generate Input Signal
t = 1:200;
xa = cos(2*pi*(1200/fs)*t) + cos(2*pi*(6000/fs)*t) + cos(2*pi*(8000/fs)*t);
% Filter input signal with designed Butterworth filter
y = filter(b,a,xa);
figure;
subplot(2,1,1)
plot(xa)
title('Un-Filtered Signal')
subplot(2,1,2)
plot(y)
title('Filtered Signal')