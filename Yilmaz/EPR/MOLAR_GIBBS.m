function out = MOLAR_GIBBS(VF,DP,X)

%% GIBBS
% Calculate Gibbs energy of a phase per mole of sites given:
% VF = [phi1 phi2]',
% DP = [DP1 DP2]',
% X = scalar
% out = Delta G of that phase per mole of sites.

out = VF/DP(1)*log(VF)+(1-VF)/DP(2)*log(1-VF) + X*VF*(1-VF);

end