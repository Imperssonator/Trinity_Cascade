function out = SINGLE_GIBBS(VF,DP,X)
%% Gibbs 1-phase 2-component
% VF = 2x1, DP = 2x1, X = scalar

out = [sum((VF./DP).*log(VF)); X*VF(1)*VF(2)];
out = [out; sum(out)];

end