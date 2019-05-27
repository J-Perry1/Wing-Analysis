function [Cl,Cd] = forces(circ,cp,delstarl,thetal,delstaru,thetau)
%Calculates the lift and drag coefficients for the aerofoil

thetate = 0.5*(thetal(length(thetal))+thetau(length(thetau))); %use average?
uete = sqrt(1-cp(length(cp))); %from equation on p27;
Hte = 0.5*(delstarl(length(delstarl))/thetal(length(thetal)) + delstaru(length(delstaru))/thetau(length(thetau)));
thetainf = thetate*(uete)^(0.5*(Hte+5));

Cd = 2*thetainf;
Cl = -2*circ;
end

