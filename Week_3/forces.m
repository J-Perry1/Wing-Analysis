function [cl,cd] = forces(circ,cp,delstarl,thetal,delstaru,thetau)
%Calculates the lift and drag coefficients for the aerofoil

thetate = thetal(length(theatal)); %use upper or lower or average?
uete = sqrt(1-cp(length(cp)); %from equation on p27;
Hte = delstarl(length(delstarl))/thetal(length(thetal));
thetainf = thetate*(uete)^(0.5*(Hte+5));

cd = 2*thetainf;
cl = -2*circ;
end

