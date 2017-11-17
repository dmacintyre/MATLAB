close all;
clear all;
clc;

n = -50 : 50;

x = repmat([1 -1],1,10);

%for i = 1:50
i = 1;
while 1
    L = i;
    top = sin((pi*n)/L);
    bottom = (pi*n)/L;

    sinc = top./bottom;

    sinc(n==0) = 1;
    clear fig;
    subplot(4,1,1)
    stem(n,sinc)
    title(L)
    x_expand = zeros(1,length(x)*L);
    x_expand(1:L:end) = x;
    y = conv(x_expand,sinc,'same');
    subplot(4,1,2)
    stem(y)
    title('Sinc Interpolation')
    linear = 1-abs(n)./L;
    linear(abs(n)>L) = 0;
    subplot(4,1,3);
    stem(n,linear);
    title(L)
    ylin = conv(linear,x_expand, 'same');
    subplot(4,1,4);
    stem(ylin);
    title('Linear Interpolation')
    pause(1) 
    if i == 50
        i = 1;
    else
        i = i + 1;
    end
end