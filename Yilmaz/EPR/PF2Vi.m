function Vi = PF2Vi(PF,VFO)
%% Phase fraction to Volumes
% convert comps x phases PF matrix to comps x phases Vi matrix which is the
% volume of each component in each phase
% constant system volume of 1 is assumed here so that we don't fuck EVERYTHING
% up

[m,n] = size(PF);
Vi = zeros(m,n);
for i = 1:n
    Vi(:,i) = VFO.*PF(:,i);
end
end