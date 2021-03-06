clear all
close all
global Re ue0 duedx
n = 101;
x = linspace(0, 1, n);
Re = 10^5;
ue0 = 1;
duedx = -0.3803;  %0.3803/0.3804 is the one for Re = 10^5
for i = 1:length(x)
    
    ue(i) = ue0 + duedx * x(i);
    
end
laminar = true;
i = 1;
f = 0;
He(1) = 1.57258; %corresponding to blasius boundary layer
itr = 0; %turbulent reattachment 
its = 0; %turbulent separation
int = 0; %natural transition
ils = 0; %laminar separation
H(1) = 0;
theta(1) = 0;
while i <= (length(x) - 1) && laminar
    i = i + 1;
    xa = x(i-1);
    xb = x(i);
    ua = ue(i-1);
    ub = ue(i);
    f = f + ueintbit(x(i-1), ue(i-1), x(i), ue(i));
    theta2 = (0.45/Re) * (ue(i))^(-6) * f;
    theta(i) = sqrt(theta2);
    %Transition calculation and check
    Rethet = Re * ue(i) * theta(i);
    m = -Re * (theta(i))^2 * duedx;
    H(i) = thwaites_lookup(m);
    He(i) = laminar_He(H(i));
    if log(Rethet) >= (18.4*He(i) - 21.74)
        laminar = false;
        int = i;

    elseif m >= 0.09
        
        laminar = false;
        ils = i; %laminar separation has occurred 
    end

end

if ils ~= 0
    He(i) = 1.51509;%set to lamianr separation value
end
deltae(i) = He(i) * theta(i);
%Setting values at the start of the panel
%As i is now at the start of the panel
while its == 0 && i <= length(x) - 1
    i = i + 1;
    %Setting initial conditions from previous iteration or before loop
    %Works as we say i = i + 1 just before the loop
    ue0 = ue(i-1);
    thick0(1) = theta(i-1);
    thick0(2) = deltae(i-1);
    [delx thicklist] = ode45(@thickdash, [0, x(i) - x(i-1)], thick0);
    
    theta(i) = thicklist(length(thicklist(:,1)), 1);
    deltae(i) = thicklist(length(thicklist(:,2)), 2);
    
    He(i) = deltae(i)/theta(i);
    %Test for reattachment
    
    if ils ~= 0 && He(i) > 1.58 && itr == 0
        itr = i;
    end
    %test for turbulent separation
    if its == 0 && He(i) < 1.46
        its = i;
    end
    
    
end
%Value of H at the turbulent separation from p.18
H = 2.803;

%calculate theta values for turbulently separated boundary layer
while i <= length(x) - 1 
    i = i + 1;
    He(i) = NaN;
    theta(i) = theta(i-1) * (ue(i-1)/ue(i))^(H+2); 
end

figure(1)
size(x)
size(He)
plot(x, He)

disp(['turbulent = ' num2str(int)])
disp(['laminar separation = ' num2str(ils)])
disp(['turbulent reaatachment = ' num2str(itr)])
disp(['turbulent separation = ' num2str(its)])