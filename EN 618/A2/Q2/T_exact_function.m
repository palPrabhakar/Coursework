function [T] = T_exact_function(x,t,alpha,sigma,x0,u)

A = -(x-x0-(u*t))^2 ;
B = sigma+(4*alpha*t) ;
C = sqrt (sigma/B) ;

T = C*(exp(A/B)) ;