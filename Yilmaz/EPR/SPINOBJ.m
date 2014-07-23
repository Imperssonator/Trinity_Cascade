function out = SPINOBJ(PF,VFO,DP)

PF = [PF,1-PF];
VF = PF2VF(PF,VFO);
obj = zeros(2,1);
G1 = Gij(VF(:,1));
G2 = Gij(VF(:,2));

% objective function: G23*G33=(G23)^2
obj(1) = (G23(VF(:,1),DP,G1)*G33(VF(:,1),DP,G1)-G23(VF(:,1),DP,G1)^2)^2;
obj(2) = (G23(VF(:,2),DP,G2)*G33(VF(:,2),DP,G2)-G23(VF(:,2),DP,G2)^2)^2;
% obj(3) = (1-PF(1,1)-PF(1,2))^2;
% obj(4) = (1-PF(2,1)-PF(2,2))^2;
% obj(5) = (1-PF(3,1)-PF(3,2))^2;

out = sum(obj);

end

