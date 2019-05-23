function [int,ils,itr,its,delstar,theta] = bl_solv(x,cp)
%UNTITLED22 Summary of this function goes here
%   Detailed explanation goes here
global Re ue0 duedx
n = length(x);

   
laminar = true;
i = 1;
f = 0;
He(1) = 1.57258; %corresponding to blasius boundary layer
itr = 0; %turbulent reattachment 
its = 0; %turbulent separation
int = 0; %natural transition
ils = 0; %laminar separation
H(1) = 0;

%initialisation
cp_init = 1;
ue(1) = 0;
    %i = i + 1;
    xa = 0;
    xb = x(1);
    ua = 0;
    ub = sqrt(1-cp(1));
    f = f + ueintbit(xa, ua, xb, ub);
    theta2 = (0.45/Re) * (ub)^(-6) * f;
    theta(1) = sqrt(theta2);
    
    %Transition calculation and check not necessary? 
 
while i <= (length(x) - 1) && laminar
    i = i + 1;
    ue(i) = sqrt(1-cp(i)); %give velocity from pressure coeffient 
    duedx = (ue(i)-ue(i-1))/(x(i)-x(i-1)); %velocity gradient 
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
    delstar(i) = H(i)*theta(i);
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
    ue(i) = sqrt(1-cp(i)); %give velocity from pressure coeffient 
    duedx = (ue(i)-ue(i-1))/(x(i)-x(i-1)); %velocity gradient 
    thick0(1) = theta(i-1);
    thick0(2) = deltae(i-1);
    [delx thicklist] = ode45(@thickdash, [0, x(i) - x(i-1)], thick0);
    
    theta(i) = thicklist(length(thicklist(:,1)), 1);
    deltae(i) = thicklist(length(thicklist(:,2)), 2);
    
    He(i) = deltae(i)/theta(i);
    %Test for reattachment
    
    %find H based on the formula on p18
    if He >= 1.46
    H = (11*He+15)/(48*He-59);
    else H = 2.803;
    end
    
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
    ue(i) = sqrt(1-cp(i)); %give velocity from pressure coeffient
    He(i) = NaN;
    delstar(i) = H*theta(i);
    theta(i) = theta(i-1) * (ue(i-1)/ue(i))^(H+2); 
end

%{
figure(1)
size(x)
size(He)
plot(x, He)

disp(['turbulent = ' num2str(int)])
disp(['laminar separation = ' num2str(ils)])
disp(['turbulent reaatachment = ' num2str(itr)])
disp(['turbulent separation = ' num2str(its)])
%}
end

