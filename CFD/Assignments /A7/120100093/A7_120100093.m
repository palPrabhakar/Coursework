clear all 
clc

h=[ 0.1, 0.05, 0.025, 0.0125, 0.00625];
x_initial=0;
x_final=5;
y_initial=0;
n=(x_final-x_initial)./h;

%exact solution 
[y1,x1]=exact(x_initial,n(1),h(1));
[y2,x2]=exact(x_initial,n(2),h(2));
[y3,x3]=exact(x_initial,n(3),h(3));
[y4,x4]=exact(x_initial,n(4),h(4));
[y5,x5]=exact(x_initial,n(5),h(5));

%RK2
[yrk21,trk2(1)]=RK2(y_initial,x1,n(1),h(1));
[yrk22,trk2(2)]=RK2(y_initial,x2,n(2),h(2));
[yrk23,trk2(3)]=RK2(y_initial,x3,n(3),h(3));
[yrk24,trk2(4)]=RK2(y_initial,x4,n(4),h(4));
[yrk25,trk2(5)]=RK2(y_initial,x5,n(5),h(5));

%RK3
[yrk31,trk3(1)]=RK3(y_initial,x1,n(1),h(1));
[yrk32,trk3(2)]=RK3(y_initial,x2,n(2),h(2));
[yrk33,trk3(3)]=RK3(y_initial,x3,n(3),h(3));
[yrk34,trk3(4)]=RK3(y_initial,x4,n(4),h(4));
[yrk35,trk3(5)]=RK3(y_initial,x5,n(5),h(5));

%RK4
[yrk41,trk4(1)]=RK4(y_initial,x1,n(1),h(1));
[yrk42,trk4(2)]=RK4(y_initial,x2,n(2),h(2));
[yrk43,trk4(3)]=RK4(y_initial,x3,n(3),h(3));
[yrk44,trk4(4)]=RK4(y_initial,x4,n(4),h(4));
[yrk45,trk4(5)]=RK4(y_initial,x5,n(5),h(5));

%AB4
[yab41,tab4(1)]=AB4(y_initial,x1,n(1),h(1)); 
[yab42,tab4(2)]=AB4(y_initial,x2,n(2),h(2));
[yab43,tab4(3)]=AB4(y_initial,x3,n(3),h(3));
[yab44,tab4(4)]=AB4(y_initial,x4,n(4),h(4));
[yab45,tab4(5)]=AB4(y_initial,x5,n(5),h(5));


%error calculation for all methods and all values of h at x=1,2,3,4,5 
err1=err(y1,yrk21,yrk31,yrk41,yab41,h(1),n(1));
err2=err(y2,yrk22,yrk32,yrk42,yab42,h(2),n(2));
err3=err(y3,yrk23,yrk33,yrk43,yab43,h(3),n(3));
err4=err(y4,yrk24,yrk34,yrk44,yab44,h(4),n(4));
err5=err(y5,yrk25,yrk35,yrk45,yab45,h(5),n(5));

%error using richardson extrapolation 
erich1=ric_err(y1,yrk21,yrk22,yrk31,yrk32,yrk41,yrk42,yab41,yab42,h(1),n(1),n(2));
erich2=ric_err(y2,yrk22,yrk23,yrk32,yrk33,yrk42,yrk43,yab42,yab43,h(2),n(2),n(3));
erich3=ric_err(y3,yrk23,yrk44,yrk33,yrk34,yrk43,yrk44,yab43,yab44,h(3),n(3),n(4));
erich4=ric_err(y4,yrk24,yrk25,yrk34,yrk35,yrk44,yrk45,yab44,yab45,h(4),n(4),n(5));

%Error in RK3 and RK4 method 
Cerr(:,:)=[err1(:,2)-err1(:,3), err2(:,2)-err2(:,3), err3(:,2)-err3(:,3), err4(:,2)-err4(:,3), err5(:,2)-err5(:,3)];

