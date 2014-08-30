function out = GIBBS2bin(PF,VFO,DP,X)
%% Gibbs 2-phase 2-component
% PF = 2x2, VFO = 2x1, DP = 2x1, X = scalar

N = PF2N(PF,VFO,DP);
Vi = PF2Vi(PF,VFO);
VF = PF2VF(PF,VFO);

out = sum(N(:,1).*log(VF(:,1))) + X*Vi(1,1)*VF(2,1)...
    + sum(N(:,2).*log(VF(:,2))) + X*Vi(1,2)*VF(2,2);

end