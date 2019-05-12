clear all
close all
np = 100;
nx = 51;
ny = 41;
xmin = -2.5;
xmax = 2.5;
ymin = -2;
ymax = 2;
theta = (0:np)*2*pi/np;
xs = cos(theta);
ys = sin(theta);
gammac = -2*sin(theta);

for i = 1:nx
    
    for j = 1:ny
        
        %Defining the grid points
        X(i,j) = xmin  + (i-1) * (xmax-xmin)/(nx-1);
        Y(i,j) = ymin  + (j-1) * (ymax-ymin)/(ny-1);
        
        psi(i,j) = Y(i,j);
        
        for k = 1:np
            
            [infa, infb] = panelinf(xs(k), ys(k), xs(k+1), ys(k+1), X(i,j), Y(i,j));
            tempor = infa * gammac(k) + infb * gammac(k+1);
            psi(i,j) = psi(i,j) + tempor;
        end
    end
    
end
    
c = -1.75:0.25:1.75;

figure(1)
contour(X,Y, psi, c)
colorbar

hold on
plot(xs,ys)
hold off