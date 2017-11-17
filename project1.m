% DSP Project 1
% Donald MacIntyre - djm4912

clear;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 60;
n = 1:M;
k = 1:4;
f1 = 0.2*pi;
p1 = 0;
f2 = 0.4*pi;
p2 = -pi/2;
f3 = 0.8*pi;
p3 = pi/5;

PART1 = 1;
PART2 = 1;

gain = 1/3.0955;

ck = 0.95*exp(j*(0.15*pi+0.02*pi*k));

h1_zeros = 0.98*exp(j*[0.8*pi -0.8*pi]'); 
h1_poles = 0.8*exp(j*[0.4*pi -0.4*pi]');

h2_zeros = [1/ck(1) 1/ck(1) 1/ck(1)' 1/ck(1)' ...
            1/ck(2) 1/ck(2) 1/ck(2)' 1/ck(2)' ...
            1/ck(3) 1/ck(3) 1/ck(3)' 1/ck(3)' ...
            1/ck(4) 1/ck(4) 1/ck(4)' 1/ck(4)']';
        
h2_poles = [ck(1) ck(1) ck(1)' ck(1)' ...
            ck(2) ck(2) ck(2)' ck(2)' ...
            ck(3) ck(3) ck(3)' ck(3)' ...
            ck(4) ck(4) ck(4)' ck(4)']';

if PART1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Pole Zero Diagram for Filter 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
zero = vertcat(h1_zeros, h2_zeros);
pole = vertcat(h1_poles, h2_poles);

figure;z = zplane(zero, pole);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Phase 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

[b,a] = zp2tf(zero,pole,gain);
[H,w] = freqz(b,a);
theta = angle(H);
unwrapped_theta = unwrap(angle(H));
t = linspace(0,pi,length(H));
t = t./pi;
figure;plot(t,radtodeg((theta)))
title('Wrapped Phase')
xlabel('Normalized Frequency');
ylabel('Phase in degrees');

figure;plot(t,radtodeg((unwrapped_theta)))
title('Unwrapped Phase')
xlabel('Normalized Frequency');
ylabel('Phase in degrees');

%figure;freqz(b,a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Group Delay 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
figure;grpdelay(b,a);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnitude of Frequency Response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
figure;plot(t,20*log(abs(H)));
title('Magnitude of Frequency Response');
xlabel('Normalized Frequency');
ylabel('Magnitude [dB]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Input Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = 0.54-0.46*cos((2*pi*n)/M);
 
x1 = w.*cos(f1*n+p1);
x2 = w.*cos(f2*n+p2);
x3 = w.*cos(f3*n+p3);

x = zeros(1,300);
x(1:M) = x3;
x(M+1:2*M) = x1;
x(2*M+1:3*M) = x2;
figure;plot(x)
title('Input Signal x[n]');
xlabel('Sample number');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot DTFT of x[n]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h1 = freqz(x);
figure;plot(t,20*log(abs(h1)));
title('DTFT of x[n]');
xlabel('Normalized Frequency');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = filter(b,a,x);
figure;plot(y)
title('Output Signal y[n]');
xlabel('Sample number');

end;

if PART2

N = 15;    
fs = 8e3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define filter H(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aa = [1 -6.6814 21.8365 -45.7544 68.2005 -76.1361 65.4106 -43.8106 22.9357 -9.3222 2.8911 -0.6625 0.1059 -0.0106 0.0005];
g = 1.91486834935659e-7;
i = [0:N-1];
b = g*[factorial(N-1)./(factorial(i).*factorial(N-1-i))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine filter stability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filt_stable = isstable(b,aa);
figure;pzmap(b,aa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine Magnitude Response of the filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h2,w2] = freqz(b,aa);

t = linspace(0,pi,length(h2));
t = t./pi;
t = t*4000;
mag = 20*log(abs(h2));
figure;plot(t,mag);
title('Filter Magnitude');
xlabel('Frequency');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the 3dB Point of the Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_mag = max(mag);
logic_three_db_down = mag<=max_mag-3;
three_db_freq = 4000*(find(logic_three_db_down,1)/length(h2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the input signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n2 = 1:100;
x2 = 1 + 2*cos(2*pi*(500/fs)*n2)+4*cos(4*pi*(10^3/fs)*n2);
figure;plot(n2,x2)
title('Input Signal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the output signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = filter(b,aa,x2);
figure;plot(y)
title('Output Signal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate ratio of input to output power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x2_squared_sum = x2*x2';
input_pwr = x2_squared_sum / length(x2);

y_squared_sum = y*y';
output_pwr = y_squared_sum / length(y);

pwr_ratio = output_pwr / input_pwr


end