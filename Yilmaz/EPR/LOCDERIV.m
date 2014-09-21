function [GRAD,DIR,E0] = LOCDERIV(PF,VFO,DP,X,ds)
%% LOCGRAD
% VFO is 2x1 [overall vol. fracs of each component]
% PF is 2x2 phase fractions
% DP is "degrees of polymerization".. 2x1
% GRAD is the vector of differential energies obtained by moving a small
% amount of each component from Phase I to Phase II.
% DIR gives the amount moved in each direction to obtain these dE's.
% Usually equal to a default step size but close to boundaries it must be
% limited to avoid INf or NaN values
% If the value is close to the boundary, the gradient will take a half step
% towards the boundary but the shell function should recognize this and
% move on

DIR = LSP(PF,ds);
possible_moves = 6;
GRAD = zeros(possible_moves,1);

E0 = GIBBS2bin(PF,VFO,DP,X);
%disp(E0)
for i = 1:possible_moves
    if any(DIR(:,i))
        GRAD(i) = GIBBS2bin([PF(:,1)+DIR(:,i),1-(PF(:,1)+DIR(:,i))],VFO,DP,X)-E0;
    else
        GRAD(i) = 1000;
    end
end

end

function DIR = LSP(PF,ds)
%% Lumpy Space Princess
% Cuts out the moves that would send you off the fuckin' MAP

DIR = [1 0 -1 0 1 -1;...
       0 1 0 -1 -1 1].*ds;

if not(any((1-PF(:,1))<ds)) && not(any(PF(:,1)<ds))
    return
end

for i = 1:2
    if PF(i,1)<ds
        badmoves = find(DIR(i,:)<0);
        DIR(:,badmoves) = zeros(2,length(badmoves));
    elseif (1-PF(i,1))<ds
        badmoves = find(DIR(i,:)>0);
        DIR(:,badmoves) = zeros(2,length(badmoves));
    end
end

end
        

