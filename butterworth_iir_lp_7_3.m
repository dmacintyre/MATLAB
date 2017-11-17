% Example 7.3 Of Text
% Butterworth Bilinear

clear all;
close;
clc;

% Define Filter Paramaters
deltap = 1-0.89125;
deltas = 0.17783;
wp = 0.2*pi;
ws = 0.3*pi;

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
figure;plot(norm_w,angle(h));

figure;grpdelay(b,a)