function [infa, infb] = refpaninf(del, X, Yin)
%This function returns influence coefficients for the given X, Yin and
%delta coefficients of a vortex sheet

%This part calculates Io

if abs(Yin) < 1e-19
    Y = 1e-19;
else
    Y = Yin;
end


Io = -(1/(4*pi)) * (X * log(X^2 + Y^2) - (X - del) * log((X-del)^2 + Y^2)...
    - 2*del + 2*Y*(atan(X/Y) - atan((X-del)/Y)));

%This part calculates I1
I1 = (1/(8*pi)) * ( (X^2 + Y^2)*log(X^2 + Y^2) - ...
    ((X-del)^2 + Y^2)*log((X-del)^2 + Y^2) - 2*X*del + del^2);


%Calculating infa, infb

tempora = (1 - (X/del)) * Io - (I1/del);

temporb = (X/del)*Io + (I1/del);

infa = tempora;
infb = temporb;
end

