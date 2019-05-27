clear all
close all
caseref = 'naca0012';
angles = -16:0.5:16;

for i = 1:length(angles)
    fname = ['Data/' caseref '_' num2str(angles(i)) '.mat'];
    load(fname)
    cd(i) = Cd;
end

plot(angles, cd)