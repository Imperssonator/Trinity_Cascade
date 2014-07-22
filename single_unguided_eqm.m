function [phase1, phase2, stability] = single_unguided_eqm(comps,DOP,Temp)

global DOPs T
DOPs = DOP; T = Temp;
balls = 1E5; %a million balls per bag
[B1,B2] = fill_bags(comps,balls) %this will create two phases (bags) with n balls specified by the compositions of each ntest
[B1e,B2e,stability] = equilibrate_bags(B1,B2) %run the equilibration algorithm to determine instability, metastability, and/eor the binodal points
phase1 = balls_to_vols(B1e);
phase2 = balls_to_vols(B2e);
end

function [B1e,B2e,stability] = equilibrate_bags(B1,B2)
%% Equilibrate Bags

% equilibrate bags takes two 1x3 vectors of integer #s of balls in two bags,
% and moves balls between those bags until the overall energy of the system
% is minimized according to a ?G function
% [B1e, B2e] are both 1x3 integer ball count vectors, stability could be:
% 0: stable - B1,B2 are at equilibrium
% 1: metastable - B1,B2 are at an energy minima but a lower energy eqm
% exists
% 2: unstable - B1,B2 spontaneously separate to an eqm energy minima

% first find if the point is at an energy minima. If not, gradient search
% to the binodal points:
% find the direction of steepest energetic descent, allowing only direction
% which are possible.
% if it is a minima, return the bags
% if not, move balls to a lower energy state and repeat

% if it is an energy minima, perform simulated annealing?

start = Gibbs_calc(B1,B2,0,0); %make the 1x7 vector of both bags and their starting energy

%find initial gradient to determine stability
[iGRAD,best_move] = get_grad(B1,B2);
disp(iGRAD)
if best_move ~=1
    Eq_bags = iterative_equilibrator([iGRAD(1,:);zeros(6,7)+1E6]); %this is a 1x7 2-bag compositions and energy vector
    B1e = Eq_bags(1:3);
    B2e = Eq_bags(4:6);
    stability = 2;
else
    energy_map(iGRAD)
    [B1e,B2e,stability] = global_min(B1,B2);
end
end

function out = energy_map(iGRAD)

E = [0 0 0];
Ei = iGRAD(1,7);
E = E+Ei;

for i = 1:3
    Energy = Ei;
    iter = 1;
    oldpt = iGRAD(1,:);
    while all(oldpt-1) && iter<10000
        iter = iter+1;
        newpt = Gibbs_calc(oldpt(1:3),oldpt(4:6),i,1);
        E(iter,i) = newpt(7);
        Energy = newpt(7);
        oldpt = newpt;
    end
end


figure
subplot(3,1,1)
plot(find(E(:,1)),E(E(:,1)~=0,1))
subplot(3,1,2)
plot(find(E(:,2)),E(E(:,2)~=0,2))
subplot(3,1,3)
plot(find(E(:,3)),E(E(:,3)~=0,3))
end

function [B1e,B2e,stability] = nucleation(iGRAD)

%% Nucleation
% input: 7x7 [GRAD] of the following format:
% [<balls of A in bag 1>... B1B B1C B2A B2B B2C <?G of this composition>]
% first row is the original test composition
% leave all 1E6's if a direction is impossible

% output: 1x7 top row of [GRAD] when the best move is no move
% iteration strategy: check if a previous iteration chose a best move and
% first try to repeat that move. If it is a good move, make it and iterate on.
% If that move is no longer a good move, go back and check all six possible
% ball moves to choose a new direction

iter = 0;
best_move = 0; % can take any value from 1 to 7

while best_move ~= 1
    %iter=iter+1;
    if best_move
        [GRAD,best_move]=make_best_move(GRAD,best_move);
        if best_move == 1 % i.e. if no move is the best move,
            [GRAD,best_move]=make_new_move(GRAD);
            if best_move == 1 % if best move is STILL no move,
                out = GRAD(1,:);
            end
        end
    else
        [GRAD,best_move]=make_new_move(GRAD);
        if best_move == 1 % if no best move and no move at all,
            out = GRAD(1,:);
        end
    end
end


end

function [B1e,B2e,stability] = global_min(B1,B2)
%% Global Min
% This runs simulated annealing on the bags because it has been determined
% that they are at a minima

global DOPs T
start_mols = B1+B2;
if all(start_mols)
    f = @(x)Gibbs_Energy(x,start_mols,DOPs,T);
    %options = saoptimset('TolFun',1e-8);
    %[eqm_mfs,Energy] = fmincon(f,[0.5 0.5 0.5],[],[],[],[],[0 0 0],[1 1 1],[],optimoptions('fmincon','Algorithm','interior-point'));
    [eqm_mfs,Energy] = simulannealbnd(f,[0.5 0.5 0.5],zeros(1,3)+1E-6,zeros(1,3)+1-1E-6,saoptimset('TolFun',1E-10));
    %disp(Energy)
else
    disp(['component ' num2str(find(start_mols==0)) 'is empty'])
end

B1e = round(start_mols.*eqm_mfs);
B2e = start_mols - B1e;

if eqm_mfs == [0.5 0.5 0.5];
    stability = 2;
else
    stability = 1;
end

end

function out = iterative_equilibrator(GRAD)
%% Iterative Equilibrator

% input: 7x7 [GRAD] of the following format:
% [<balls of A in bag 1>... B1B B1C B2A B2B B2C <?G of this composition>]
% first row is the original test composition
% leave all 1E6's if a direction is impossible

% output: 1x7 top row of [GRAD] when the best move is no move
% iteration strategy: check if a previous iteration chose a best move and
% first try to repeat that move. If it is a good move, make it and iterate on.
% If that move is no longer a good move, go back and check all six possible
% ball moves to choose a new direction

iter = 0;
best_move = 0; % can take any value from 1 to 7

while best_move ~= 1
    %iter=iter+1;
    if best_move
        [GRAD,best_move]=make_best_move(GRAD,best_move);
        if best_move == 1 % i.e. if no move is the best move,
            [GRAD,best_move]=make_new_move(GRAD);
            if best_move == 1 % if best move is STILL no move,
                out = GRAD(1,:);
            end
        end
    else
        [GRAD,best_move]=make_new_move(GRAD);
        if best_move == 1 % if no best move and no move at all,
            out = GRAD(1,:);
        end
    end
end

end

function [GRAD,best_move] = get_grad(B1,B2)

GRAD = [Gibbs_calc(B1,B2,0,0); zeros(6,7)];

for i = 1:3
    if B1(i)>1 %if there are any balls of i in B1,
        GRAD(2*i,:) = Gibbs_calc(B1,B2,i,1); %find the energy of moving an i from bag 1
    else
        GRAD(2*i,:) = zeros(1,7)+1E6;
    end
    if B2(i)>1 %if there are any i in B2,
        GRAD(2*i+1,:) = Gibbs_calc(B1,B2,i,2); %find energy of moving i from bag
    else
        GRAD(2*i+1,:) = zeros(1,7)+1E6;
    end
end

[lowest_E, best_move] = min(GRAD(:,7)); %this will output the lowest energy value and its index

end

function [GRAD,best_move] = make_best_move(oldGRAD,move)
%% Make best move
%take a 6x7 [oldGRAD] and make 'move' (scalar, 2-7), then decide if that move
%is the best move, i.e. best_move will be scalar, 1 OR 'move'

MOVES = [0 0;...
         1 1;...
         1 2;...
         2 1;...
         2 2;...
         3 1;...
         3 2];
comp = MOVES(move,1); %the ball we move
from = MOVES(move,2); %the bag we move it from
GRAD = [oldGRAD(1,:); zeros(1,7)];

if oldGRAD(1,comp+3*(from-1))>1 % if there are >1 balls of 'comp' in bag 'from'
    GRAD(2,:) = Gibbs_calc(GRAD(1,1:3),GRAD(1,4:6),comp,from);
else
    GRAD(2,:) = zeros(1,7)+1E6;
end

[lowest_E, best_move] = min(GRAD(:,7));
GRAD = [GRAD(best_move,:); zeros(6,7)+1E6];
    
end

function [GRAD,best_move] = make_new_move(oldGRAD)

%% Make new move
% make new move computes the energy of moving any of the three components
% to or from each bag, 6 moves in total, returning [GRAD] with the first
% row as the new best, its previous index as the best move, and large
% positive numbers to fill the rest.
% 2: ball 1 from bag 1, 3: ball 1 from bag 2
% 4: ball 2 from bag 1, 5: ball 2 from bag 2
% 6: ball 3 from bag 1, 7: ball 3 from bag 2
GRAD = [oldGRAD(1,:); zeros(1,7)];

for i = 1:3
    if oldGRAD(1,i)>1 %if there are any balls of i in B1,
        GRAD(2*i,:) = Gibbs_calc(oldGRAD(1,1:3),oldGRAD(1,4:6),i,1); %find the energy of moving an i from bag 1
    else
        GRAD(2*i,:) = zeros(1,7)+1E6;
    end
    if GRAD(1,i+3)>1 %if there are any i in B2,
        GRAD(2*i+1,:) = Gibbs_calc(GRAD(1,1:3),GRAD(1,4:6),i,2); %find energy of moving i from bag
    else
        GRAD(2*i+1,:) = zeros(1,7)+1E6;
    end
end

[lowest_E, best_move] = min(GRAD(:,7)); %this will output the lowest energy value and its index
GRAD = [GRAD(best_move,:); zeros(6,7)+1E6];

end

function out = Gibbs_calc(B1,B2,comp,from)

%% Gibbs Calc

% Gibbs calc takes two 1x3 vectors that are the bags, changes the ball
% counts to reflect which *component (scalar, 1 2 or 3) is being moved *from what bag (scalar, 1 or 2),
% calculates the energy of the new conformation, and returns this as a 1x7
% vector:
% [B1A B1B B1C B2A B2B B2C ?G]
% pass a '0' to 'from' to calculate energy without switching balls around

global T

R = 8.314;
M1 = B1; %start with the new bags the same as given
M2 = B2;

%the following section moves a ball according to which component it is and
%which bag it's taken from

if from == 1
    M1(comp) = B1(comp) - 1;
    M2(comp) = B2(comp) + 1;
elseif from == 2
    M1(comp) = B1(comp) + 1;
    M2(comp) = B2(comp) - 1;
end

%now we convert total ball counts to mol fractions, then volume fractions
VF1 = balls_to_vols(M1);
VF2 = balls_to_vols(M2);

Energy1 = sum(M1.*log(VF1)) + g12()*M1(1)*VF1(2) + g13()*M1(1)*VF1(3) + g23()*M1(2)*VF1(3);
Energy2 = sum(M2.*log(VF2)) + g12()*M2(1)*VF2(2) + g13()*M2(1)*VF2(3) + g23()*M2(2)*VF2(3);
Energy = R*T*(Energy1+Energy2);
%disp(Energy)

% actual formula: ?G/RT = n1*log(vf1) + n2*log(vf2) + n3*log(vf3) +
% g12(u2)*n1*vf2 + g13(vf3)*n1*vf3 + g23(vf3)*n2vf3
% in this case we are using the ball counts as number of moles. Dividing
% by total number of moles would yield molar gibbs but we're just going to add
% it all up because I'm lazy

out = [M1, M2, Energy];

end

function VF = balls_to_vols(B)

global DOPs

denom = 0;
for i = 1:3
    denom = denom + B(i)*DOPs(i);
end

for i = 1:2
    VF(i) = B(i)*DOPs(i)/denom;
end

VF(3) = 1 - VF(1) - VF(2);

end

function out = g12()
%out = Xij('CHCl3','P3HT',295);
%out = 0.99; %P3HT/CHCl3
%out = 0.1;
out = -0.3; %water/acetone
end

function out = g13()
%out = Xij('CHCl3','PS',295);
%out = 1.5; %P3HT/"hexane"
%out = 0.16;
out = 1.0; %water/cellulose acetate
end

function out = g23()
%out = 0.05; %CHCl3/"hexane"
%out = 0.004;
out = 0.2; %acetone/cellulose acetate
end

function comps = vol_to_mol(vf)
%takes a 1x3 vector and outputs another 1x3 vector

global DOPs

denom = 0; %first calculate the denominator of the formula
for i = 1:3
    denom = denom + vf(i)/DOPs(i);
end

for i = 1:2
    comps(i) = (vf(i)/DOPs(i))/denom;
end

comps(3) = 1 - comps(1) - comps(2);

end

function [B1,B2] = fill_bags(VF,balls)

B1 = zeros(1,3); %each bag will have a certain number of balls

comps = vol_to_mol(VF); %the starting compositions are volume fractions and must be converted to mol fractions here to fill the bags
B1(1) = ceil(balls*comps(1)); %multiply mol fractions by total balls in each bag
B1(2) = floor(balls*comps(2));
B1(3) = balls - B1(1) - B1(2); %fill polymers first, then solvent

B2 = B1; %both bags start with the same number of balls

end