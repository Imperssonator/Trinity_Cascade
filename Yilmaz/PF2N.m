function N = PF2N(PF,VFO,DP)
%% Phase fraction to Moles
% convert comps x phases PF matrix to comps x phases N matrix which is the
% number of moles in each phase
% constant system volume is assumed here so that we don't fuck EVERYTHING
% up

[m,n] = size(PF);
N = zeros(m,n);
for i = 1:n
    N(:,i) = VFO.*PF(:,i)./DP;
end
end