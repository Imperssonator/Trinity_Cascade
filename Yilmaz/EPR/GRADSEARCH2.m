function [PF,E,stab] = GRADSEARCH2(PFi,VFO,DP,pert,initstep,conTol,species,T)

VFi = PF2VF(PFi,VFO);
X = Chi(species,T)+pert;
[GRi,DIRi,Ei] = LOCDERIV(PFi,VFO,DP,X,1E-6);
%disp(GRi)

if not(any(GRi<0)) % if no initial directions decrease energy
    PF = PFi;
    E = Ei;
    stab = 1; % give up the gradient search
else
    [PF,E] = downgrad(PFi,VFO,DP,pert,Ei,initstep,conTol,species,T);
    stab = 0;
end

end

function [PF,E,stab] = downgrad(PFi,VFO,DP,pert,Ei,step,conTol,species,T)
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
X = Chi(species,T)+pert;
[GR,DIR,E] = LOCDERIV(PFi,VFO,DP,X,step); % initialize the gradient with the new step size
PF = PFi;
if not(any(GR<0))
    stab = 1; % if the gradient is positive in all directions, terminate
end

hist_size = 40;
move_thresh = 5;
MHIST = zeros(2,hist_size);
EHIST = zeros(1,hist_size);
EHIST(end) = E;
iter = 0;
freq = 4;
super_move_min = 5; %multiple of previous move's energy difference that the super move must beat to be used

while stab == 0
    iter = iter+1;
    made_move = 0;
    if ( freq * round(double(iter)/freq) == iter ) %every 'freq' moves
        %disp('testing turbo')
        if max(abs(sum(MHIST')))>move_thresh*step % if the move history sum is greater than a threshold in a particular direction
            [is_good_move,dE] = try_move(PF,VFO,DP,X,sum(MHIST')',EHIST,super_move_min); %check to see where repeating the whole move history would go
            if is_good_move
                move = sum(MHIST')';
                %disp('turbo!!')
                %disp(move)
                made_move = 1;
            end
        end
    end
    if made_move == 0
        [move,dE] = pick_best_move(GR,DIR); % picks the move vector in the direciton of greatest energy decrease
    end
%     disp('_________')
%     disp(PF)
%     disp(move)
%     disp(dE)
    PF = make_move(PF,move); % applies the move vector to the current phase fraction
    VF = PF2VF(PF,VFO);
    X = Chi(species,T)+pert;
    MHIST = [MHIST(:,2:hist_size),move];
    EHIST = [EHIST(2:hist_size),E+dE];
    [GR,DIR,E] = LOCDERIV(PF,VFO,DP,X,step); % checks the gradient at the new point
    %disp(GR)
    if not(any(GR<0))
        %disp(move)
        stab = 1;
    end
end

[PF,E,stab] = downgrad(PF,VFO,DP,pert,E,step/2,conTol,species,T);

end

function [is_good_move,dE] = try_move(PF,VFO,DP,X,move,EHIST,super_move_min)

PFnew = zeros(2,2);
PFnew(:,1) = PF(:,1)+move;
PFnew(:,2) = 1-PFnew(:,1);


if all([PFnew(:,1);PFnew(:,2)]>0) && all([PFnew(:,1);PFnew(:,2)]<1)
    dE = GIBBS2bin(PFnew,VFO,DP,X)-EHIST(end);
    %disp('dE for turbo:')
    %disp(dE)
    if dE/super_move_min > (EHIST(end)-EHIST(end-1))
        is_good_move = 0;
        
    else
        is_good_move = 1;
    end
else
    is_good_move = 0;
    dE = 1E6;
    %disp('turbo overshot')
end

end

function PF = make_move(PFi,move)

PF(:,1) = PFi(:,1)+move;
PF(:,2) = 1-PF(:,1);
%disp(PF)

end

function [move,best_dE] = pick_best_move(GR,DIR)
% out is a move: 2x1

[best_dE,best_dir] = min(GR);
move = DIR(:,best_dir);

end