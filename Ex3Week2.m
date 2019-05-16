clear all
close all
%This script plots theta for a flat plate using the solution described in
%the handout (Ex1 week 2) as well as compares it with the Blasius method
ReL = linspace(10e3, 10e7, 100000); 
ueL = [0.5];
n = 101;
theta = zeros(length(ReL), length(ueL), n);
control=1;
%Initialising first point
for i = 1:length(ReL)
    for j = 1:length(ueL)
        x = linspace(0, 1, n);
        ue = linspace(1, ueL(j), n);
        laminar = true;
        k = 1;
        f = 0;
        int = 0;
        ils = 0;
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
                int = k+1;
                if control == 1
                    Re_transsep = ReL(i);
                    control = 0;
                end
                
            
            elseif m >= 0.09
                
                laminar = false;
                ils = k+1;
            end
            
            k = k + 1;
        end
        end
        
        if int ~= 0
           % disp(['Natural transition at ' num2str(x(int))...
           %     ' with Rethet ' num2str(Rethet)])
            
        end
        
        if ils ~= 0
            %    disp([' Separation at ' num2str(x(ils))...
            %    ' with Rethet ' num2str(Rethet)])
        end
           
        
end

disp(['Surpass ReL = ' num2str(Re_transsep)])

