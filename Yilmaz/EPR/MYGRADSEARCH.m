function [PF,E,stab] = MYGRADSEARCH(PFi,VFO,DP,pert,initstep,conTol)

VFi = PF2VF(PFi,VFO);
G1 = Gij(VFi(:,1))+pert;
G2 = Gij(VFi(:,2))+pert;
[GRi,DIRi,Ei] = LOCGRAD(PFi,VFO,DP,G1,G2,initstep);

if not(any(GRi<0))
    PF = PFi;
    E = Ei;
    stab = 1;
else
    [PF,E] = downgrad(PFi,VFO,DP,pert,Ei,initstep,conTol);
    stab = 0;
end

end

function [PF,E,stab] = downgrad(PFi,VFO,DP,pert,Ei,step,conTol)
%% DOWNGRAD
% Ride the gradient downhill by a certain step size until you hit a minima.
% Then decrease the step size and repeat until the step size is smaller
% than conTol
% a 'move' is a 3x1 vector that would be added to PF1 and subtracted from
% PF2

stab = 0; % if this has been called, it is unstable
if step < conTol % if we're taking small enough steps, abort
    PF = PFi;
    E = Ei;
    stab = 1;
    return
end

VFi = PF2VF(PFi,VFO);
G1 = Gij(VFi(:,1))+pert;
G2 = Gij(VFi(:,2))+pert;
[GR,DIR,E] = LOCGRAD(PFi,VFO,DP,G1,G2,step); % initialize the gradient with the new step size
PF = PFi;
if not(any(GR<0))
    stab = 1; % if the gradient is positive in all directions, terminate
end

% hist_size = 25;
% move_thresh = 15;
% MHIST = zeros(3,hist_size);
% iter = 0;

while stab == 0
%     iter = iter+1;
%     if max(sum(MHIST'))>move_thresh
%         
    move = pick_best_move(GR,DIR); % picks the move vector in the direciton of greatest energy decrease
    %disp(move)
    PF = make_move(PF,move); % applies the move vector to the current phase fraction
    VF = PF2VF(PF,VFO);
    G1 = Gij(VF(:,1))+pert;
    G2 = Gij(VF(:,2))+pert;
    [GR,DIR,E] = LOCGRAD(PF,VFO,DP,G1,G2,step); % checks the gradient at the new point
    if not(any(GR<0))
        stab = 1;
    end
end

[PF,E,stab] = downgrad(PF,VFO,DP,pert,E,step/2,conTol);

end

function PF = make_move(PFi,move)

PF(:,1) = PFi(:,1)+move;
PF(:,2) = 1-PF(:,1);
%disp(PF)

end

function out = pick_best_move(GR,DIR)
% out is a move: 3x1

[best_dE,best_dir] = min(GR);
out = DIR(:,best_dir);

end