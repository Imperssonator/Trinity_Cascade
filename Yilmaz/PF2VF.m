function VF = PF2VF(PF,VFO)
%% Phase Fraction to Volume Fraction
% Convert comps x phases PF matrix and overall vol. frac VFO vector to
% comps x phases VF matrix

[m,n] = size(PF);
VF = zeros(m,n);
for i = 1:n
    VF(:,i) = VFO.*PF(:,i)/dot(VFO,PF(:,i));
end

end