function [Xs,Ys] = unit_vect_panel(xa, ya, xb, yb, x, y)
%This function calculates influence coefficients a and b from a general vortex sheet using the previous...
%function refpaninf


%Getting tangential unit vector

vt = [xb-xa, yb-ya];
ut = vt/norm(vt);

%Getting normal unit vector (by deducting orthogonal vector from dot
%product)

vn = [-(yb-ya), xb - xa];
un = vn/norm(vn);

%Define a position vector and new coordinates
r = [x-xa, y-ya];
Xs = dot(r, ut);
Ys = dot(r,un);
end