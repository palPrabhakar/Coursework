function [yd] = der(x,h);

yd = zeros(3,2);

% First Order Scheme 

yd(1,1) = (f(x+h)-f(x-h))/(2*h); %Central 
yd(2,1) = (f(x+h) - f(x))/h; %Forward
yd(3,1) = (f(x)-f(x-h))/h; % Backward

% Second  Order Scheme 

yd(1,2) = (-f(x+2*h)+8*f(x+h)-8*f(x-h)+f(x-2*h))/(12*h); %Central 
yd(2,2) = (4*f(x+h)-f(x+2*h)-3*f(x))/(2*h); %Forward
yd(3,2) = (3*f(x)-4*f(x-h)+f(x-2*h))/(2*h); % Backward