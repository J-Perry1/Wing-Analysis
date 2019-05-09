close all 
clear all
%Scripts makes arrays of xm, ym and psi and later plots contours of psi on
%a described field
xmin = -2.5;
xmax = 2.5;
nx = 51;
ymin = -2;
ymax = 2;
ny = 41; 
xc = 0.75;
yc = 0.5;
gamma = 3;

for i = 1:nx
    
    for j = 1:ny
        
        xm(i,j) = xmin  + (i-1) * (xmax-xmin)/(nx-1);
        ym(i,j) = ymin  + (j-1) * (ymax-ymin)/(ny-1);
        psi(i,j) = psipv(xc, yc, gamma, xm(i,j), ym(i,j));
        
    end
    
end

c = -0.4:0.2:1.2;
figure
contour(xm, ym, psi, c)