clear all
close all
%This script plots theta for a flat plate using the solution described in
%the handout (Ex1 week 2) as well as compares it with the Blasius method
ReL = 10^5; 
x = linspace(0, 1, 101);
ue = linspace(1, 1, 101);
%Initialising first point
f = 0;
theta(1) = 0;
for i = 1:(length(x) - 1)
    
    xa = x(i);
    xb = x(i+1);
    ua = ue(i);
    ub = ue(i+1);
    f = f + ueintbit(x(i), ue(i), x(i+1), ue(i+1));
    theta2 = (0.45/ReL) * (ue(i))^(-6) * f;
    theta(i+1) = sqrt(theta2);
end


%Calculating Blasius solution

thetab = (0.664/(ReL^0.5)) * x.^(0.5); %x. multiplies by every element in the array

figure(1)
plot(x, theta, '-.')
hold on
plot(x, thetab)
xlabel('$x/L$','Interpreter','Latex')
ylabel('$\theta/L$','Interpreter','Latex')
legend('Thwaites'' Solution', 'Blasius'' Solution')
legend('show')