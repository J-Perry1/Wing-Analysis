clear all
close all
global Re ue0 duedx
n = 101;
x = linspace(0, 1, n);
Re = 10^5;
ue0 = 1;
gradients = linspace(0, -0.8, 1000); 
%0.3803/0.3804 is the one for Re = 10^5
control = 1;

for j = 1:length(gradients)
    ue0 = 1;
    duedx = gradients(j);
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
            ils = i;
        end

    end

    if ils ~= 0
        He(i) = 1.51509;
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

        if its == 0 && He(i) < 1.46
            its = i;
        end


    end
    %Value of H at the turbulent separation from p.18
    H = 2.803;

    while i <= length(x) - 1
        i = i + 1;
        He(i) = NaN;
        theta(i) = theta(i-1) * (ue(i-1)/ue(i))^(H+2); 
    end
    if control == 1 && its == n
        control = j;
    end
end

disp(gradients(control))
