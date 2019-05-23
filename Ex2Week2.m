clear all
close all
control = 0;
%This script plots theta for a flat plate using the solution described in
%the handout (Ex1 week 2) as well as compares it with the Blasius method
ReL = (10^6) * [1, 10, 100]; 
ueL = [0.8, 1, 1.2];
n = 101;
theta = zeros(length(ReL), length(ueL), n);
p = 1;
%Initialising first point
for i = 1:length(ReL)
    for j = 1:length(ueL)
        x = linspace(0, 1, n);
        ue = linspace(1, ueL(j), n);
        laminar = true;
        k = 1;
        f = 0;
        %disp([j ueL(j)])
        while k <= (length(x) - 1) && laminar 
            xa = x(k);
            xb = x(k+1);
            ua = ue(k);
            ub = ue(k+1);
            f = f + ueintbit(x(k), ue(k), x(k+1), ue(k+1));
            theta2 = (0.45/ReL(i)) * (ue(k+1))^(-6) * f;
            theta(i,j,k+1) = sqrt(theta2);
            
            %Transition calculation and check
            Rethet = ReL(i) * ue(k+1) * theta(i,j, k+1);
            m = -ReL(i) * (theta(i,j,k+1))^2 * (ue(length(ue)) - ue(1));
            H = thwaites_lookup(m);
            He = laminar_He(H);
            if log(Rethet) >= (18.4*He - 21.74)
                laminar = false;
                disp(ReL(i)/(10^6));
                disp(ueL(j));
                disp([x(k+1), Rethet/1000]);
                tab_x(p) = x(k+1);
                tab_Rethet(p) = Rethet/1000;
                p = p+1;
                %disp(k)
            end
            
            k = k + 1;
            end
            
        end
end
for i = 1:n
    theta11(i) = theta(1,1,i);
    theta12(i) = theta(1,2,i);
    theta13(i) = theta(1,3,i);
    theta21(i) = theta(2,1,i);
    theta22(i) = theta(2,2,i);
    theta23(i) = theta(2,3,i);
    theta31(i) = theta(3,1,i);
    theta32(i) = theta(3,2,i);
    theta33(i) = theta(3,3,i);
end
plot(x, theta21)
hold on
plot(x, theta22)
hold on
plot(x, theta23)
legend('0.8', '1', '1.2')


tab_x = transpose(tab_x);
tab_Rethet = transpose(tab_Rethet);

T = table(tab_x,tab_Rethet);

%Calculating Blasius solution

%thetab = (0.664/(ReL^0.5)) * x.^(0.5);

%figure(1)
%plot(x, theta)
%hold on
%plot(x, thetab)
%xlabel('x/L')
%ylabel('theta/L')
%legend('Thwaites solution', 'Blasius solution')
%legend('show')