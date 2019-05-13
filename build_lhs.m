function lhsmat = build_lhs(xs,ys)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
np = length(xs) -1;
psip = zeros(np,np+1);
lhsmat = zeros(np+1,np+1);
for i = 1:np
    for j = 1:np+1
        X(i) = (xs(i)+xs(i+1))/2; %vortex influence is at midpoint? 
        Y(i) = (ys(i)+ys(i+1))/2;
        [infa(j),infb(j)] = panelinf(xs(i), ys(i), xs(i+1), ys(i+1), X(i), Y(i));
        if j==1
            psip(i,j) = infa(j);
        elseif j == np+1
            psip(i,j) = infb(j-1);
        else 
            psip(i,j) = infa(j) + infb(j-1);
        end 
        
        if i == j && i < np %diagonal matrix
            lhsmat(i,j) = lhsmat(i,j) + psip(i+1,j) - psip(i,j);
        end 
    end
end

end

