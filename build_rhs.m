function rhsvec = build_rhs(xs,ys,alpha)
%build_rhs function to build b matrix of equation 9

np = length(xs)-1;
%initialise the matrices to ensure correct size 
psifs = zeros(np,1);
rhsvec = zeros(np+1,1);
for i = 1:np+1 
    psifs(i) = ys(i)*cos(alpha)-xs(i)*sin(alpha);
end

for i = 1:np-1 
    rhsvec(i) = psifs(i) - psifs(i+1);
end
end
