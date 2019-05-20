
close all
clear all
global Re ue0 duedx;

U = 1;
Re_table = [10^6, 10^7, 10^8];

for i = 1:length(Re_table)
    Re = Re_table(i);
    ue0 = U;
    duedx = -0.6
    control = 0;
    x0 = 0.01;
    thick0(1) = 0.037*x0*(Re*x0)^(-1/5);
    thick0(2) = 1.8*thick0(1);

    [delx, thickhist] = ode45(@thickdash,[0 0.99],thick0);
    
    for j = 1:length(delx)
        x(j,i) = x0 + delx(j);
    end
    
    for j = 1:length(delx)
        deltae(j,i) = thickhist(j,2);
        theta(j,i) = thickhist(j,1);
    end
    
    for j = 1:length(x)
        He(j, i) = deltae(j, i)/theta(j, i);
            if He(j, i) < 1.46 && control == 0 && x(j, i) > 0.05
                xsep(i) = x(j, i);
                control = 1
            end
    end
    
end
x(x==0) = NaN;
He(He==0) = NaN;
theta(theta==0) = NaN;
deltae(deltae==0) = NaN;
figure(1)
plot(x(:,1), He(:,1), 'DisplayName', 'Re = 10^6')
hold on
plot(x(:,2), He(:,2), 'DisplayName', 'Re = 10^7')
hold on
plot(x(:,3), He(:,3), 'DisplayName', 'Re = 10^8')
xlabel('x/L')
ylabel('He')
ylim([0, 2.5])
yline(1.46, 'DisplayName', 'Separation Point')
title('Plot of He against x/L for grad = -0.6')
legend('show')

figure(2)
%plot(x(:,1), theta(:,1), 'DisplayName', 'Theta for Re = 10^6')
hold on
plot(x(:,2), theta(:,2), 'DisplayName', 'Theta for Re = 10^7 and grad = -0.6')
hold on
%plot(x(:,3), theta(:,3), 'DisplayName', 'Theta for Re = 10^8')
hold on
plot(x(:,2), deltae(:,2), 'DisplayName', 'DeltaE for Re = 10^7 and grad = -0.6')
% plot(x(:,1), deltae(:,1), 'DisplayName', 'DeltaE for grad = -0.3')
% hold on
% plot(x(:,2), deltae(:,2), 'DisplayName', 'DeltaE for grad = -0.6')
% hold on
% plot(x(:,3), deltae(:,3), 'DisplayName', 'DeltaE for grad = -0.9')
xlabel('x/L')
legend('show')

disp(xsep)
