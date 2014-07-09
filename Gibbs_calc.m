function out = Gibbs_calc(X,moles,DOPs,T)

%% Gibbs Calc
% 

VF = mols_to_vols([moles, X],DOPs);
VF1 = VF(1:3)
VF2 = VF(4:6)
M1 = moles.*X
M2 = moles.*(1-X)

R = 8.314;

out = R*T*(sum(M1.*log(VF1)) + g12()*M1(1)*VF1(2) + g13()*M1(1)*VF1(3) + g23()*M1(2)*VF1(3) ...
    + sum(M2.*log(VF2)) + g12()*M2(1)*VF2(2) + g13()*M2(1)*VF2(3) + g23()*M2(2)*VF2(3));

% actual formula: ?G/RT = n1*log(vf1) + n2*log(vf2) + n3*log(vf3) +
% g12(u2)*n1*vf2 + g13(vf3)*n1*vf3 + g23(vf3)*n2vf3

end

function out = mols_to_vols(M,DOPs)

%% mols to vols
% takes a 1x6 [total moles in system, molar fraction of each comp. in phase
% 1] and outputs another 1x6, [VFs phase 1, VFs phase 2]

M1 = M(1:3).*M(4:6);
M2 = M(1:3) - M1;

out(1:3) = M2V(M1,DOPs);
out(4:6) = M2V(M2,DOPs);

end

function out = M2V(M,DOPs)
% just an inner function of mols to vols that needs to be run twice,
% converts 1x3 vector of moles to 1x3 vector of VFs

denom = 0;
for i = 1:3
    denom = denom+M(i)*DOPs(i);
end
for i = 1:2
    out(i) = M(i)*DOPs(i)/denom;
end
out(3) = 1 - out(1) - out(2);
end

function out = g12()
%out = Xij('CHCl3','P3HT',295);
out = 0.99; %CHCl3/P3HT
%out = 0.4;
end

function out = g13()
%out = Xij('CHCl3','PS',295);
out = 0.39; % CHCl3/PS
%out = 0.4;
end

function out = g23()
out = 0.48; %PS/P3HT
%out = 0.004;
end