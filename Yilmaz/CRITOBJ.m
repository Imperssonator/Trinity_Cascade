function out = CRITOBJ(VF,DP)

obj = zeros(2,1);
G = Gij(VF);
g22 = G22(VF,DP,G);
g23 = G23(VF,DP,G);

obj(1) = (1-DP(1)/DP(2)*(VF(1)/VF(2))^2-2*g22/g23*(1-g22/g23)-1-DP(1)/DP(3)*(VF(1)/VF(3))^2*(g22/g23)^3)^2;
obj(2) = (1-sum(VF))^2;
out = sum(obj);
end