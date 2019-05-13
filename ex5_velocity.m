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
alpha = 0.8; %angle of attack

A = build_lhs(xs,ys);
b = build_rhs(xs,ys,alpha);
gam = A\b; %solve the matrix equation
theta_plot = theta/pi; %normalise theta for the plot

circ = 0;
%calculate the total circulation by summing values of gamma x r x dtheta 
for i = 1:length(gam)
    circ = circ + gam(i)*2*pi/np;
end
disp(circ)

figure 
plot(theta_plot,gam)

axis([0 2 -2.5 2.5])
