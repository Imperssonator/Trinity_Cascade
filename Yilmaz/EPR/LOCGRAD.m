function [GRAD,DIR,E0] = LOCGRAD(PF,VFO,DP,G1,G2,ds)
%% LOCGRAD
% VFO is 3x1 [overall vol. fracs of each component]
% PF is 3x2 phase fractions
% DP is "degrees of polymerization".. 3x1
% GRAD is the vector of differential energies obtained by moving a small
% amount of each component from Phase I to Phase II.
% DIR gives the amount moved in each direction to obtain these dE's.
% Usually equal to a default step size but close to boundaries it must be
% limited to avoid INf or NaN values
% If the value is close to the boundary, the gradient will take a half step
% towards the boundary but the shell function should recognize this and
% move on

default_step = ds;
DIR = largest_step_possible(PF,default_step);
GRAD = zeros(6,1);

E0 = GIBBS2(PF,VFO,DP,G1,G2);
for i = 1:6
    GRAD(i) = GIBBS2([PF(:,1)+DIR(:,i),PF(:,2)-DIR(:,i)],VFO,DP,G1,G2)-E0;
end

end

function DIR = largest_step_possible(PF,ds)

DIR = zeros(3,6);
step_frac = 2;
for i = 1:3
    if PF(i,1)<ds
        DIR(i,i+3) = -PF(i,1)/step_frac; DIR(i,i) = ds;
    elseif 1-PF(i,1)<ds
        DIR(i,i) = (1-PF(i,1))/step_frac; DIR(i,i+3) = -ds;
    else
        DIR(i,i) = ds; DIR(i,i+3) = -ds;
    end
end
end
        

