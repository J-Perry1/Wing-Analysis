function psixy = psipv(xc, yc, gamma, x, y)
% Functions returns a value of the stream function of the vortex centred at
% (xc, yc) of strength gamma at point (x,y) See ex.1 handout

r2 = (x-xc)^2 + (y-yc)^2;
psi = - (gamma/(4*pi)) * log(r2);
psixy = psi;
end

