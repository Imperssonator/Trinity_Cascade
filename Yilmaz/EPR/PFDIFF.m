function out = PFDIFF(PF)
%% sum squared phase fraction difference

PFD = zeros(6,1);
PFD(1) = PF(1,1)-PF(2,1);
PFD(2) = PF(1,1)-PF(3,1);
PFD(3) = PF(2,1)-PF(3,1);
PFD(4) = PF(1,2)-PF(2,2);
PFD(5) = PF(1,2)-PF(3,2);
PFD(6) = PF(2,2)-PF(3,2);

out = sum(PFD.^2);

end