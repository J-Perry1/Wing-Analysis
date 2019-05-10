close all 
clear all

xmin = 0;
xmax = 5;
nx = 51;
ymin = 0;
ymax = 4;
ny = 41;
del = 1.5;
nv = 100;
xa = 4.1;
ya = 1.3;
xb = 2.2;
yb = 2.9;
for i = 1:nx
    
    for j = 1:ny
        
        %Defining the grid points
        X(i,j) = xmin  + (i-1) * (xmax-xmin)/(nx-1);
        Y(i,j) = ymin  + (j-1) * (ymax-ymin)/(ny-1);
        %Calculating infa and infb using the external function
        [infa(i,j), infb(i,j)] = panelinf(del, xa, ya, xb, yb, X(i,j), Y(i,j));
        
    end
    
end
%c = -0.15:0.05:0.15;

figure(1)
contourf(X, Y,infa)
%title('infa formula - general sheet')
colorbar

figure(2)
contourf(X, Y,infb)
%title('infb formula - general sheet')
colorbar
