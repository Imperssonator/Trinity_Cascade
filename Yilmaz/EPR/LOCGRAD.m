function [GRAD,step] = LOCGRAD(PF,VFO,DP)
%% LOCGRAD
% VFO is 3x1 [overall vol. fracs of each component]
% PF is 3x2 phase fractions
% DP is "degrees of polymerization".. 3x1
% GRAD is the vector of differential energies obtained by moving a small
% amount of each component from Phase I to Phase II.
% the small amount is given by 'step', which is a default value or whatever
% distance brings the system only a fraction of a step closer to the edge
% of phase space

default_step = 0.001;
G = Gij(VFO);
DIR = eye(3).*default_step;
GRAD = zeros(3,1);

E0 = GIBBS2(PF,VFO,DP,G,G);
for i = 1:3
    GRAD(i) = GIBBS2([PF(:,1)+DIR(:,i),PF(:,2)-DIR(:,i)],VFO,DP,G,G)-E0;
end

end