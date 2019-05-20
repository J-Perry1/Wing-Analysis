
close all
clear all
global Re ue0 duedx;

Re = 5*10^6;
ue0 = 1;
duedx = 0;

x0 = 0.01;
thick0(1) = 0.037*x0*(Re*x0)^(-1/5);
thick0(2) = 1.8*thick0(1);

[delx, thickhist] = ode45(@thickdash,[0 0.99],thick0);

deltae = thickhist(:,2);
theta = thickhist(:,1);

%power law fits
%x = linspace(x0,0.99,length(theta);
x = x0 + delx;
for i = 1:length(x)
    theta7(i) = 0.037*x(i)*(Re*x(i))^(-1/5); %1/7th power law
    theta9(i) = 0.023*x(i)*(Re*x(i))^(-1/6); %1/9th power law 
end

figure 
plot(x,theta,'DisplayName','DE Solution')
hold on
plot(x,theta7,'DisplayName','1/7th Power Law')
plot(x,theta9,'DisplayName','1/9th Power Law')
legend
