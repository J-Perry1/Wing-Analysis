
%This script plots influence coefficients a and b for the vortices from the
%expressions derived in the notes as well as approximates it using the
%finite vortex sheet small vortex contributions
xmin = -2.5;
xmax = 2.5;
nx = 51;
ymin = -2;
ymax = 2;
ny = 41;
del = 1.5;
nv = 100;
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
            if abs(Y(i,j)) < 1e-19
                Yin = 1e-19;
            else
                Yin = Y(i,j);
            end
            tempor_psi = psipv(L, 0, gam, X(i,j), Yin);
            psi_a(i,j) = psi_a(i,j) + tempor_psi;
            
        end
        
        psi_b(i,j) = 0;
        for k = 0:nv
            L = k*del/nv;
            gamx = k * (1/nv); %vulnerable to error
            gam = gamx * (del/nv);
            if abs(Y(i,j)) < 1e-19
                Yin = 1e-19;
            else
                Yin = Y(i,j);
            end
            tempor_psi = psipv(L, 0, gam, X(i,j), Yin);
            psi_b(i,j) = psi_b(i,j) + tempor_psi;
        end
        
    end
    
end
c = -0.15:0.05:0.15;

figure(1)
contour(X, Y,infa,c)
%title('infa formula')
xlabel('X')
ylabel('Y')
colorbar
h = colorbar;
ylabel(h, 'Streamfunction Value')

figure(2)
contour(X, Y,psi_a, c)
%title('infa discrete vortex approximation')
xlabel('X')
ylabel('Y')
colorbar
h = colorbar;
ylabel(h, 'Streamfunction Value')

figure(3)
contour(X, Y,infb, c)
%title('infb formula')
xlabel('X')
ylabel('Y')
colorbar
h = colorbar;
ylabel(h, 'Streamfunction Value')

figure(4)
contour(X, Y,psi_b, c);
%title('infb discrete vortex approximation')
colorbar
h = colorbar;
ylabel(h, 'Streamfunction Value')