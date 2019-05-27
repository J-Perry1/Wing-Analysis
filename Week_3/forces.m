function [Cl,Cd] = forces(circ,cp,delstarl,thetal,delstaru,thetau)
%Calculates the lift and drag coefficients for the aerofoil

thetate = thetal(end)+thetau(end); %use average?
uete = sqrt(1-cp(end)); %from equation on p27;
Hte = (delstarl(end) + delstaru(end))/thetate;
thetainf = thetate*(uete)^(0.5*(Hte+5));

Cd = 2*thetainf;
Cl = -2*circ;
end

