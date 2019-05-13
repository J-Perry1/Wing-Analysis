function lhsmat = build_lhs(xs,ys)
%build_lhs function to build A matrix for equation 9
%   Detailed explanation goes here
np = length(xs) -1;
psip = zeros(np,np+1);
lhsmat = zeros(np+1,np+1);
for i = 1:np
    for j = 1:np
        %X(i) = (xs(i)+xs(i+1))/2; %vortex influence is at midpoint? 
        %Y(i) = (ys(i)+ys(i+1))/2;
        [infa(j),infb(j)] = panelinf(xs(j), ys(j), xs(j+1), ys(j+1), xs(i), ys(i));
        if j==1
            psip(i,j) = infa(j);
        else 
            psip(i,j) = infa(j) + infb(j-1);
        end 
    end
end
j = np+1;
[infa(j),infb(j)] = panelinf(xs(j-1), ys(j-1), xs(j), ys(j), xs(i), ys(i));
psip(i,j) = infb(j-1);
disp(psip)
for i = 1:np
    for j = 1:np+1
        if i == j && i < np %diagonal matrix
            lhsmat(i,j) = lhsmat(i,j) + psip(i+1,j) - psip(i,j);
        end 
    end
end
A(np+1,np+1) = 1;
A(1,1) = 1;
end

