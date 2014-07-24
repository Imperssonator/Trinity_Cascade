function out = GIBBS(N,VF,Vi,G)

%% GIBBS
% Calculate Gibbs energy of a phase given:
% N = [n1 n2 n3], VF = [phi1 phi2 phi3],
% Gij = 3x3 G(1,2) = g12, G(1,3) = g13, G(2,3) = g23
% out = Delta G of that phase. Total, not per mole.

out = sum(N.*log(VF)) + G(1,2)*Vi(1)*VF(2) + G(1,3)*Vi(1)*VF(3) + G(2,3)*Vi(2)*VF(3);

end