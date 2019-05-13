clear all
close all
np = 100;
nx = 51;
ny = 41;
xmin = -2.5;
xmax = 2.5;
ymin = -2;
ymax = 2;
theta = (0:np)*2*pi/np;
xs = cos(theta);
ys = sin(theta);
gammac = -2*sin(theta);
alpha = 0;

A = build_lhs(xs,ys);
b = build_rhs(xs,ys,alpha);
gam = A\b;
theta_plot = theta/pi;
%c = -1.75:0.25:1.75;

figure 
plot(theta_plot,gam)
axis([0 2 -2.5 2.5])
