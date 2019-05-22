clear all
close all
global Re ue0 duedx
n = 101;
x = linspace(0, 1, n);
Retab = [10^4, 10^5, 10^6];
ue0 = 1; 
%0.3803/0.3804 is the one for Re = 10^5
control = 0;

for j = 1:length(Retab)
    ue0 = 1;
    duedx = -0.25;
    Re = Retab(j);
    for i = 1:length(x)

        ue(i) = ue0 + duedx * x(i);

    end
    laminar = true;
    i = 1;
    f = 0;
    He(j,1) = 1.57258;
    itr(j) = 0;
    theta_tr(j) = 0;
    He_tr(j) = 0;
    its(j) = 0;
    theta_ts(j) = 0;
    He_ts(j) = 0;
    int(j) = 0;
    theta_nt(j) = 0;
    He_nt(j) = 0;
    ils(j) = 0;
    theta_ls(j) = 0;
    He_ls(j) = 0;
    H(1) = 0;
    theta(j,1) = 0;
    while i <= (length(x) - 1) && laminar
        i = i + 1;
        xa = x(i-1);
        xb = x(i);
        ua = ue(i-1);
        ub = ue(i);
        f = f + ueintbit(x(i-1), ue(i-1), x(i), ue(i));
        theta2 = (0.45/Re) * (ue(i))^(-6) * f;
        theta(j,i) = sqrt(theta2);

         %Transition calculation and check
        Rethet = Re * ue(i) * theta(j,i);
        m = -Re * (theta(j,i))^2 * duedx;
        H(i) = thwaites_lookup(m);
        He(j,i) = laminar_He(H(i));
        if log(Rethet) >= (18.4*He(j,i) - 21.74)
            laminar = false;
            int(j) = i;
            theta_nt(j) = theta(j,i);
            He_nt(j) = He(j,i);

        elseif m >= 0.09

            laminar = false;
            ils(j) = i;
            theta_ls(j) = theta(j,i);
            He_ls(j) = He(j,i);
        end

    end

    if ils(j) ~= 0
        He(j,i) = 1.51509;
    end
    deltae(i) = He(j,i) * theta(j,i);
    %Setting values at the start of the panel
    %As i is now at the start of the panel
    while its(j) == 0 && i <= length(x) - 1
        i = i + 1;
        %Setting initial conditions from previous iteration or before loop
        %Works as we say i = i + 1 just before the loop
        ue0 = ue(i-1);
        thick0(1) = theta(j,i-1);
        thick0(2) = deltae(i-1);
        [delx thicklist] = ode45(@thickdash, [0, x(i) - x(i-1)], thick0);

        theta(j,i) = thicklist(length(thicklist(:,1)), 1);
        deltae(i) = thicklist(length(thicklist(:,2)), 2);

        He(j,i) = deltae(i)/theta(j,i);
        %Test for reattachment

        if ils(j) ~= 0 && He(j,i) > 1.58 && itr(j) == 0
            itr(j) = i;
            theta_tr(j) = theta(j,i);
            He_tr(j) = He(j,i);
        end

        if its(j) == 0 && He(j,i) < 1.46
            its(j) = i;
            theta_ts(j) = theta(j,i);
            He_ts(j) = He(j,i);
        end


    end
    %Value of H at the turbulent separation from p.18
    H = 2.803;

    while i <= length(x) - 1
        i = i + 1;
        He(j,i) = NaN;
        theta(j,i) = theta(j,i-1) * (ue(i-1)/ue(i))^(H+2); 
    end
    

end

    figure(1)
    plot(x, theta(1,:), 'DisplayName', 'Re = 10^4')
    hold on
    plot(x, theta(2,:), 'DisplayName', 'Re = 10^5')
    hold on
    plot(x, theta(3,:), 'DisplayName', 'Re = 10^6')
    if its(1) ~= 0
        hold on
        scatter(x(its(1)), theta_ts(1), 'DisplayName', 'Turbulent separation', 'Marker', 'o', 'MarkerFaceColor', 'r')
    end
    
    if its(2) ~= 0
        hold on
        scatter(x(its(2)), theta_ts(2), 'Marker', 'o', 'MarkerFaceColor', 'r')
    end
    
    if its(3) ~= 0
        hold on
        scatter(x(its(3)), theta_ts(3), 'Marker', 'o', 'MarkerFaceColor', 'r')
    end
    
    if ils(1) ~= 0
        hold on
        scatter(x(ils(1)), theta_ls(1), 'DisplayName', 'Laminar separation', 'Marker', 'o', 'MarkerFaceColor', 'b')
    end
    
    if ils(2) ~= 0
        hold on
        h_special_1 = scatter(x(ils(2)), theta_ls(2), 'Marker', 'o', 'MarkerFaceColor', 'b')
        set( get( get( h_special_1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    end
    
    if ils(3) ~= 0
        hold on
        scatter(x(ils(3)), theta_ls(3), 'Marker', 'o', 'MarkerFaceColor', 'b')
    end
    
    if int(1) ~= 0
        hold on
        scatter(x(int(1)), theta_nt(1), 'Marker', 'o', 'MarkerFaceColor', 'c')
    end
    
    if int(2) ~= 0
        hold on
        scatter(x(int(2)), theta_nt(2),'Marker', 'o', 'MarkerFaceColor', 'c')
    end
    
    if int(3) ~= 0
        hold on
        scatter(x(int(3)), theta_nt(3), 'DisplayName', 'Natural transition', 'Marker', 'o', 'MarkerFaceColor', 'c')
    end
    
    if itr(1) ~= 0
        hold on
        scatter(x(itr(1)), theta_tr(1), 'DisplayName', 'Turbulent reattachment', 'Marker', 'o', 'MarkerFaceColor', 'm')

    end
    
    if itr(2) ~= 0
        hold on
        scatter(x(itr(2)), theta_tr(2), 'DisplayName', 'Turbulent reattachment', 'Marker', 'o', 'MarkerFaceColor', 'm')
    end
    
    if itr(3) ~= 0
        hold on
        scatter(x(itr(3)), theta_tr(3), 'DisplayName', 'Turbulent reattachment', 'Marker', 'o', 'MarkerFaceColor', 'm')
    end
    xlabel('x/L')
    ylabel('theta/L')
    ylim([0 0.022])
    legend
    
    %Here starts plot 2
    figure(2)
    plot(x, He(1,:), 'DisplayName', 'Re = 10^4')
    hold on
    plot(x, He(2,:), 'DisplayName', 'Re = 10^5')
    hold on
    plot(x, He(3,:), 'DisplayName', 'Re = 10^6')
    if its(1) ~= 0
        hold on
        scatter(x(its(1)), He_ts(1), 'DisplayName', 'Turbulent separation', 'Marker', 'o', 'MarkerFaceColor', 'r')
    end
    
    if its(2) ~= 0
        hold on
        scatter(x(its(2)), He_ts(2),  'DisplayName', 'Turbulent separation', 'Marker', 'o', 'MarkerFaceColor', 'r')
    end
    
    if its(3) ~= 0
        hold on
        scatter(x(its(3)), He_ts(3),  'DisplayName', 'Turbulent separation', 'Marker', 'o', 'MarkerFaceColor', 'r')
    end
    
    if ils(1) ~= 0
        hold on
        scatter(x(ils(1)), He_ls(1),  'DisplayName', 'Laminar separation', 'Marker', 'o', 'MarkerFaceColor', 'b')
    end
    
    if ils(2) ~= 0
        hold on
        h_special_2 = scatter(x(ils(2)), He_ls(2), 'Marker', 'o', 'MarkerFaceColor', 'b')
        set( get( get( h_special_2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    end
    
    if ils(3) ~= 0
        hold on
        scatter(x(ils(3)), He_ls(3),  'DisplayName', 'Laminar separation', 'Marker', 'o', 'MarkerFaceColor', 'b')
    end
    
    if int(1) ~= 0
        hold on
        scatter(x(int(1)), He_nt(1),  'DisplayName', 'Natural transition', 'Marker', 'o', 'MarkerFaceColor', 'c')
    end
    
    
    if int(2) ~= 0
        hold on
        scatter(x(int(2)), He_nt(2), 'DisplayName', 'Natural transition', 'Marker', 'o', 'MarkerFaceColor', 'c')
    end
    
    if int(3) ~= 0
        hold on
        scatter(x(int(3)), He_nt(3),  'DisplayName', 'Natural transition', 'Marker', 'o', 'MarkerFaceColor', 'c')
    end
    
    if itr(1) ~= 0
        hold on
        scatter(x(itr(1)), He_tr(1), 'DisplayName', 'Turbulent reattachment', 'Marker', 'o', 'MarkerFaceColor', 'm')

    end
    
    if itr(2) ~= 0
        hold on
        scatter(x(itr(2)), He_tr(2), 'DisplayName', 'Turbulent reattachment', 'Marker', 'o', 'MarkerFaceColor', 'm')
    end
    
    if itr(3) ~= 0
        hold on
        scatter(x(itr(3)), He_tr(3), 'DisplayName', 'Turbulent reattachment', 'Marker', 'o', 'MarkerFaceColor', 'm')
    end
    
    xlabel('x/L')
    ylabel('He')
    ylim([0 2])
    legend('Location', 'SouthEast')