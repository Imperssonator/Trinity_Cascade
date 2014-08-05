function [GRAD,DIR,E0] = LOCGRAD(PF,VFO,DP,G1,G2,ds)
%% LOCGRAD
% VFO is 3x1 [overall vol. fracs of each component]
% PF is 3x2 phase fractions
% DP is "degrees of polymerization".. 3x1
% GRAD is the vector of differential energies obtained by moving a small
% amount of each component from Phase I to Phase II.
% DIR gives the amount moved in each direction to obtain these dE's.
% Usually equal to a default step size but close to boundaries it must be
% limited to avoid INf or NaN values
% If the value is close to the boundary, the gradient will take a half step
% towards the boundary but the shell function should recognize this and
% move on

DIR = LSP(PF,ds);
GRAD = zeros(12,1);

E0 = GIBBS2(PF,VFO,DP,G1,G2);
%disp(E0)
for i = 1:12
    if any(DIR(:,i))
        GRAD(i) = GIBBS2([PF(:,1)+DIR(:,i),PF(:,2)-DIR(:,i)],VFO,DP,G1,G2)-E0;
    else
        GRAD(i) = 1000;
    end
end

end

function DIR = LSP(PF,ds)
%% Lumpy Space Princess
% Cuts out the moves that would send you off the fuckin' MAP

DIR = [1 0 0 -1 0 0 1 -1 0 0 1 -1;...
       0 1 0 0 -1 0 -1 1 1 -1 0 0;...
       0 0 1 0 0 -1 0 0 -1 1 -1 1].*ds;

if not(any((1-PF(:,1))<ds)) && not(any(PF(:,1)<ds))
    return
end

for i = 1:3
    if PF(i,1)<ds
        badmoves = find(DIR(i,:)<0);
        DIR(:,badmoves) = zeros(3,length(badmoves));
    elseif (1-PF(i,1))<ds
        badmoves = find(DIR(i,:)>0);
        DIR(:,badmoves) = zeros(3,length(badmoves));
    end
end

end
        

