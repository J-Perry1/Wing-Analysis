global Re ue0 duedx
n = 101;
x = linspace(0, 1, n);
Re = 10^6;
ue0 = 1;
duedx = -0.9;
for i = 1:length(x)
    
    ue(i) = ue0 + duedx * x(i);
    
end
laminar = true;
i = 1;
f = 0;
He(1) = 1.57258;
itr = 0;
its = 0;
int = 0;
ils = 0;
theta(1) = 0;
while i <= (length(x) - 1) && laminar
    xa = x(i);
    xb = x(i+1);
    ua = ue(i);
    ub = ue(i+1);
    f = f + ueintbit(x(i), ue(i), x(i+1), ue(i+1));
    theta2 = (0.45/Re) * (ue(i+1))^(-6) * f;
    theta(i+1) = sqrt(theta2);

    %Transition calculation and check
    Rethet = Re * ue(i+1) * theta(i+1);
    m = -Re * (theta(i+1))^2 * duedx;
    H = thwaites_lookup(m);
    He = laminar_He(H);
    if log(Rethet) >= (18.4*He - 21.74)
        laminar = false;
        int = i+1;

    elseif m >= 0.09

        laminar = false;
        ils = i+1;
    end

    i = i + 1;
end

if ils ~= 0
    He = 1.51509;
end

deltae(i) = He * theta(i) %As i = i + 1 at the end of the laminar loop
%Setting values at the start of the panel
%As i is now at the start of the panel
i = i + 1;
while its == 0 && i <= length(x)
    %Setting initial conditions from previous iteration or before loop
    %Works as we say i = i + 1 just before the loop
    ue0 = ue(i-1)
    thick0(1) = theta(i-1);
    thick0(2) = deltae(i-1)
    [delx thicklist] = ode45(@thickdash, [0, x(i) - x(i-1)], thick0)
    
    theta(i) = thicklist(length(thicklist(:,1)), 1)
    deltae(i) = thicklist(length(thicklist(:,2)), 2)
    
    He(i) = deltae(i)/theta(i);
    %Test for reattachment
    
    if ils ~= 0 && He(i) > 1.58 && itr == 0
        itr = i;
    end
    
    if its == 0 && He(i) < 1.46
        its = i;
    end
    
    i = i + 1;
    
end
%Value of H at the turbulent separation from p.18
H = 2.803;

while i < length(x)
    He(i) = NaN
    theta(i) = theta(i-1) * (ue(i-1)/ue(i))^(H+2) 
    
    i = i +1;
end

disp(['turbulent = ' num2str(int)])
disp(['laminar separation = ' num2str(ils)])
disp(['turbulent reaatachment = ' num2str(itr)])
disp(['turbulent separation = ' num2str(its)])
