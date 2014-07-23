function out = GIBBS2(PF,VFO,DP,G1,G2)
%% Gibbs 2-phase
% PF = 3x2, VFO = 3x1, DP = 3x1, G1, G2 = 3x3

N = PF2N(PF,VFO,DP);
VF = PF2VF(PF,VFO);

out = GIBBS(N(:,1),VF(:,1),G1)+GIBBS(N(:,2),VF(:,2),G2);

end