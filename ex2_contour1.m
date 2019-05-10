
%This script plots influence coefficients a and b for the vortices from the
%expressions derived in the notes as well as approximates it using the
%finite vortex sheet small vortex contributions

close all 
clear all

xmin = -2.5;
xmax = 2.5;
nx = 51;
ymin = -2;
ymax = 2;
ny = 41;
del = 1.5;
nv = 1000;
for i = 1:nx
    
    for j = 1:ny
        
        %Defining the grid points
        X(i,j) = xmin  + (i-1) * (xmax-xmin)/(nx-1);
        Y(i,j) = ymin  + (j-1) * (ymax-ymin)/(ny-1);
        %Calculating infa and infb using the external functions
        [infa(i,j), infb(i,j)] = refpaninf(del, X(i,j), Y(i,j));
        
        %Calculating the discrete vortex sheet contributions to get inf a
        %equivalent (gamma_a = 1, gamma_b = 0
        %zeroing psi_a
        psi_a(i,j) = 0;
        for k = 0:nv
            L = k*(del/nv);
            gamx = 1 - k * (1/nv); %vulnerable to error
            gam = gamx * (del/nv);
            tempor_psi = psipv(L, 0, gam, X(i,j), Y(i,j));
            psi_a(i,j) = psi_a(i,j) + tempor_psi;
            
        end
        
        psi_b(i,j) = 0;
        for k = 0:nv
            L = k*del/nv;
            gamx = k * (1/nv); %vulnerable to error
            gam = gamx * (del/nv);
            tempor_psi = psipv(L, 0, gam, X(i,j), Y(i,j));
            psi_b(i,j) = psi_b(i,j) + tempor_psi;
        end
        
    end
    
end
c = -0.15:0.05:0.15;

figure(1)
contourf(X, Y,infa,c)
%title('infa formula')
colorbar

figure(2)
contourf(X, Y,psi_a, c)
%title('infa discrete vortex approximation')
colorbar

figure(3)
contourf(X, Y,infb, c)
%title('infb formula')
colorbar

figure(4)
contourf(X, Y,psi_b, c);
%title('infb discrete vortex approximation')
colorbar