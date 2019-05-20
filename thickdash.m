function dthickdx = thickdash(xmx0,thick)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here
global Re ue0 duedx; %needed to evaluate differential equation

ue = ue0 + duedx*xmx0;
He = thick(2)/thick(1);

theta = thick(1);
Rethet = Re * ue * theta;

if He >= 1.46
    H = (11*He+15)/(48*He-59);
else H = 2.803;
end

 
Cf = 0.091448*((H-1)*Rethet)^(-0.232)*exp(-1.26*H);
Cdiss = 0.010019*((H-1)*Rethet)^(-1/6);

top = 0.5*Cf - ((H+1)/ue)*duedx*thick(1);
bottom = Cdiss - (3/ue)*duedx*thick(2);
f = [top,bottom];
f = transpose(f);
dthickdx = f;
end

