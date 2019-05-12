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
gammac = -2sin(theta);

for i = 1:nx
    
    for j = 1:ny
        
        %Defining the grid points
        X(i,j) = xmin  + (i-1) * (xmax-xmin)/(nx-1);
        Y(i,j) = ymin  + (j-1) * (ymax-ymin)/(ny-1);
        
        for 
        
        end
    end