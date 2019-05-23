function lhsmat = build_lhs(xs,ys)
%lhsmat builds matrix A from equation 9 using equations 5,6,7,8
np = length(xs) - 1;
%initialise matrices to ensure right size
psip = zeros(np, np+1); 
A = zeros(np+1, np+1);

for i = 1:np %iterate over rows
    
    j = 1; %special condition for first column
    [infaj , infbj] = panelinf(xs(j), ys(j), xs(j+1), ys(j+1), xs(i), ys(i));
    psip(i,j) = infaj;
    for j = 2:np %iterate over columns 
        [infaj_1, infbj_1] = panelinf(xs(j-1), ys(j-1), xs(j), ys(j), xs(i), ys(i));
        [infaj , infbj] = panelinf(xs(j), ys(j), xs(j+1), ys(j+1), xs(i), ys(i));
        psip(i,j) = infaj + infbj_1;
  
   
    end
    j = np + 1; %special condition for last column 
    [infaj_1, infbj_1] = panelinf(xs(j-1), ys(j-1), xs(j), ys(j), xs(i), ys(i));
    psip(i,j) = infbj_1;
    
end


for i = 1:np-1 %make rows of A using rows of psip via equation 6
    A(i,:) = psip(i+1,:) - psip(i,:);
end
%boundary conditions provide last two equations to solve for all unknowns
%Interpolation from code completion 1
A(np,1) = 1;
A(np, 2) = -1;
A(np,3) = 0.5;
A(np, np-1) = - 0.5;
A(np, np) = 1;
%Second equation (gam np+1 + gam 1 = 0 cause velocities on both sides are
%equal
A(np + 1, 1) = 1
A(np+1,np+1) = 1;
lhsmat = A;
end