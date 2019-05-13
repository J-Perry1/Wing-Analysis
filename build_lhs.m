function lhsmat = build_lhs(xs,ys)

np = length(xs) - 1;
psip = zeros(np, np+1);
A = zeros(np+1, np+1);

for i = 1:np
    
    j = 1;
    [infaj , infbj] = panelinf(xs(j), ys(j), xs(j+1), ys(j+1), xs(i), ys(i));
    psip(i,j) = infaj;
    for j = 2:np
        [infaj_1, infbj_1] = panelinf(xs(j-1), ys(j-1), xs(j), ys(j), xs(i), ys(i));
        [infaj , infbj] = panelinf(xs(j), ys(j), xs(j+1), ys(j+1), xs(i), ys(i));
        psip(i,j) = infaj + infbj_1;
  
   
    end
    j = np + 1;
    [infaj_1, infbj_1] = panelinf(xs(j-1), ys(j-1), xs(j), ys(j), xs(i), ys(i));
    psip(i,j) = infbj_1;
    
end


for i = 1:np-1
    A(i,:) = psip(i+1,:) - psip(i,:);
end
A(np,1) = 1;
A(np+1,np+1) = 1;
lhsmat = A;
end