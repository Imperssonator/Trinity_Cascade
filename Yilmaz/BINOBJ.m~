function out = BINOBJ(PF,VFO,DP,pert)

PF = [PF,1-PF];
VF = PF2VF(PF,VFO);
N = PF2N(PF,VFO,DP);

obj = 0;
G1 = Gij(VF(:,1))+pert;
G2 = Gij(VF(:,2))+pert;

obj = GIBBS(N(:,1),VF(:,1),G1)+GIBBS(N(:,2),VF(:,2),G2);
% obj(1) = (CHEMPOT1(VF(:,1),DP,G1)-CHEMPOT1(VF(:,2),DP,G2))^2;
% obj(2) = (CHEMPOT2(VF(:,1),DP,G1)-CHEMPOT2(VF(:,2),DP,G2))^2;
% obj(3) = (CHEMPOT3(VF(:,1),DP,G1)-CHEMPOT3(VF(:,2),DP,G2))^2;
% obj(4) = (1-PF(1,1)-PF(1,2))^2;
% obj(5) = (1-PF(2,1)-PF(2,2))^2;
% obj(6) = (1-PF(3,1)-PF(3,2))^2;

out = obj;

end