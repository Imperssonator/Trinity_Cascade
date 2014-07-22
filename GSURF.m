function out = GSURF()

pts = populate_triangle(100);
for i = 1:length(pts)
    pts(i,4) = GIBBS_E(pts(i,1:3));
end

ternsurf(pts(:,1),pts(:,2),pts(:,4))

out = pts;

end

function out = GIBBS_E(VF)

%% GIBBS
% Calculate Gibbs energy of a phase given:
% N = [n1 n2 n3], VF = [phi1 phi2 phi3],
% Gij = 3x3 G(1,2) = g12, G(1,3) = g13, G(2,3) = g23
% out = Delta G of that phase. Total, not per mole.

G = [0  0.004 0.4;...
     0  0     0.4;...
     0  0     0];
 
N = vf2N(VF); 

out = sum(N.*log(VF)) + G(1,2)*N(1)*VF(2) + G(1,3)*N(1)*VF(3) + G(2,3)*N(2)*VF(3);

end

function comps = vf2N(vf)
%takes a 1x3 vector and outputs another 1x3 vector

molar_vols = [1000 1000 1]; %molar volumes need not be absolute, just relative [non-solvent solvent polymer]

denom = 0; %first calculate the denominator of the formula
for i = 1:3
    denom = denom + vf(i)/molar_vols(i);
end

for i = 1:2
    comps(i) = (vf(i)/molar_vols(i))/denom;
end

comps(3) = 1 - comps(1) - comps(2);

end

function realout = populate_triangle(vertices)

%this takes a number of vertices and outputs an nx3 matrix of compositions
%which are volume fractions

ntests=0; %define how many starting solution concentrations we will use. vertices are how much we will segment each side of the ternary diagram
for i = 1:vertices
    ntests=ntests+i;
end

out = zeros(ntests,3);

A = linspace(.001,0.999,vertices);
B = linspace(.001,0.999,vertices);
C = linspace(.001,0.999,vertices);

i = 0; %this will count up every time we start a new line
col = 1; %this will count up every time we make a point

while col <= ntests
    for j = 1:(vertices-i)
        out(col,1) = A(j);
        out(col,2) = B(i+1);
        out(col,3) = 1-out(col,1)-out(col,2);
        col = col+1;
    end
    i = i+1;
end

row = 1;
for i = 1:length(out)
    if all(out(i,:)>0)
        realout(row,:) = out(i,:);
        row = row+1;
    end
end

end 