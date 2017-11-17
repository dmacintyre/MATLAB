% Project 2 - Q2
% Donald MacIntyre - djm4912

clear all;
close all;
clc;

% Filter Paramaters
ws = 0.25*pi;
wp = 0.375*pi;
wc = (ws+wp)/2;
deltaw = wp - ws;
delta = 0.01;
A = -20*log10(0.01);

% Determine Filter Paramaters of beta and M
beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);     % A of 40 so use range specified by 21 <= A <= 50
M_pure = (A-8) / (2.285*deltaw);
M = ceil(M_pure);
alpha = M/2;

% Generate the window with the found paramaters
n = 0:1:M;
w = besseli(0,beta*sqrt(1-((n-alpha)/alpha).^2))/besseli(0,beta);
% figure;stem(w);
% title('Kaiser Window')

% Check that designed window matches ideal Matlab window
w_matlab = kaiser(M+1,beta)';
err = max(abs(w-w_matlab));

% Create ideal highpass filter only for samples that will have a corresponding 
% non-zero window value 
nideal = -M/2:1:M/2;
hideal = (sin(pi*(nideal)))./(pi*(nideal))- (sin(wc*(nideal))./(pi*(nideal)));
% Use l'hopital's rule to calculate term when n = 0 (i.e. 0/0)
hideal(19) = 1-wc / pi;

h = hideal .*w;

%figure;stem(h)
%title('Kaiser Filter Taps')

% Plot frequency response of Filter
% [hfilt,wfilt] = freqz(h);
% wfilt = wfilt / pi;
% figure;plot(wfilt,abs(hfilt))
% title('Magnitude Response')
% xlabel('Normalized Frequency')
% ylabel('Amplitude')

% Redesign Filter to meet missed spec
M = M + 2;
alpha = M/2;
n = 0:1:M;
w = besseli(0,beta*sqrt(1-((n-alpha)/alpha).^2))/besseli(0,beta);

% Create ideal highpass filter only for samples that will have a corresponding 
% non-zero window value 
nideal = -M/2:1:M/2;
hideal = (sin(pi*(nideal)))./(pi*(nideal))- (sin(wc*(nideal))./(pi*(nideal)));
% Use l'hopital's rule to calculate term when n = 0 (i.e. 0/0)
hideal(20) = 1-wc / pi;

h = hideal .*w;

figure;stem(h)
title('Kaiser Filter Taps')

[hfilt,wfilt] = freqz(h);
wfilt = wfilt / pi;
figure;plot(wfilt,abs(hfilt))
title('Magnitude Response')
xlabel('Normalized Frequency')
ylabel('Amplitude')

