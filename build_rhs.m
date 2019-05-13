function rhsvec = build_rhs(xs,ys,alpha)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
np = length(xs)-1;
psifs = zeros(np,1);
rhsvec = zeros(np+1,1);
for i = 1:np+1
    psifs(i) = ys(i)*cos(alpha)-xs(i)*sin(alpha);
    if i < np 
        rhsvec(i) = psifs(i+1) - psifs(i);
    end
end

