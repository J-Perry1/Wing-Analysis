
close all
clear all
global Re ue0 duedx;

U = 1;

Re = 10^8;
ue0 = U;
duedx = -0.6;
%duedx = [-0.3,-0.6,-0.9];

x0 = 0.01;
thick0(1) = 0.037*x0*(Re*x0)^(-1/5);
thick0(2) = 1.8*thick0(1);

[delx, thickhist] = ode45(@thickdash,[0 0.99],thick0);

deltae = thickhist(:,2);
theta = thickhist(:,1);

x = x0 + delx;

for i = 1:length(theta)
    He(i) = deltae(i)/theta(i);
end

figure
plot(x,He)




