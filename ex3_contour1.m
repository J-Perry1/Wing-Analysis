close all 
clear all

xmin = 0;
xmax = 5;
nx = 51;
ymin = 0;
ymax = 4;
ny = 41;
nv = 500; %resolution of 100 gives artefact for infa - suggest 500
xa = 4.1;
ya = 1.3;
xb = 2.2;
yb = 2.9;
del = [xb-xa, yb-ya];
del = norm(del);

for i = 1:nx
    
    for j = 1:ny
        
        %Defining the grid points
        X(i,j) = xmin  + (i-1) * (xmax-xmin)/(nx-1);
        Y(i,j) = ymin  + (j-1) * (ymax-ymin)/(ny-1);
        %Calculating infa and infb using the external function
        [infa(i,j), infb(i,j)] = panelinf(xa, ya, xb, yb, X(i,j), Y(i,j));
        
        [Xs, Ys] = unit_vect_panel(xa, ya, xb, yb, X(i,j), Y(i,j));
        psi_a(i,j) = 0; 
        for k = 0:nv
            L = k*(del/nv);
            gamx = 1 - k * (1/nv); %vulnerable to error
            gam = gamx * (del/nv);
            if abs(Ys) < 1e-19 
                Yin = 1e-19;
            else
                Yin = Ys;
            end
            tempor_psi = psipv(L, 0, gam, Xs, Yin);
            psi_a(i,j) = psi_a(i,j) + tempor_psi;
            
        end
        
        psi_b(i,j) = 0;
        for k = 0:nv
            L = k*del/nv;
            gamx = k * (1/nv); %vulnerable to error
            gam = gamx * (del/nv);
            if abs(Ys) < 1e-19
                Yin = 1e-19;
            else
                Yin = Ys;
            end
            tempor_psi = psipv(L, 0, gam, Xs, Yin);
            psi_b(i,j) = psi_b(i,j) + tempor_psi;
        end
        
    end
    
end


%c = -0.15:0.05:0.15;

figure(1)
contour(X, Y,infa)
%title('infa formula - general sheet')
colorbar

figure(2)
contour(X, Y,infb)
%title('infb formula - general sheet')
colorbar

figure(3)
contour(X, Y,psi_a)
%title('infa discrete vortex approximation general sheet')
colorbar

figure(4)
contour(X, Y,psi_b);
%title('infb discrete vortex approximation general sheet')
colorbar
