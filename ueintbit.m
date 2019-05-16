function f = ueintbit(xa,ua,xb,ub)
%ueintbit returns a value for the integral part of Thwaites' method for
%solution of momentum integral equation
ubar = 0.5*(ua+ub);
deltau = ub-ua;
deltax = xb-xa;

f = (ubar^5+(5/6)*(ubar^3)*(deltau^2)+(1/16)*ubar*(deltau^4))*deltax;
end

