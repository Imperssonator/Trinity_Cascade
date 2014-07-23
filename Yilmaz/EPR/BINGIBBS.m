function out = BINGIBBS(VF)

DP = [1 50];
G = 2;

out = VF/DP(1)*log(VF) + (1-VF)/DP(2)*log(1-VF) + G*VF/DP(1)*(1-VF);

end