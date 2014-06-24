function [A,B,C] = evaporative_ternary()
tic
[A,B,C] = construct_binodal(50);
figure
ternplot(A,C,B,'or')
ternlabel('CHCl3','PS','P3HT')
toc
end

function [A,B,C] = construct_binodal(steps)

%% Construct Binodal
% this function takes an integer 'steps', and performs a gradient search
% for equilibrium phase compositions starting from solutions along a line that follows solvent evaporation

% A, B, and C are equal sized vectors of the compositions of each component
% to be fed to the ternplot function.

ntests = populate_triangle(steps); %create n starting solutions that follow a line from a dilute solution to a pure even mixture of two polymers

balls = 50000; %50000 balls per bag
test = 1; %count how many test points we've used
i = 1; %count how many points on the binodal we've marked
ABC = zeros(1,3); %ABC will store A, B and C as columns in a tests x 3 matrix
figure
hold on

while test<=steps
    
    [B1,B2] = fill_bags(ntests(test,:),balls); %this will create two phases (bags) with n balls specified by the compositions of each ntest
    %disp(B1)
    %disp(B2)
    
    [B1e,B2e] = equilibrate_bags(B1,B2); %run the equilibration algorithm
    %disp(B1e)
    %disp(B2e)
    VF = balls_to_vols(B1);
    
    if B1==B1e
        %disp('point was stable')
        %disp(test)
        plot(VF(1),VF(2),'ob',VF(1),VF(3),'oc')
        test = test+1; %if the point is already equilibrated, keep going
        
    else
        ABC(i,:) = balls_to_vols(B1e); %else mark two points in A, B, and C, one from each phase
        plot(VF(1),ABC(i,2),'or',VF(1),ABC(i,3),'om')
        i = i+1;
        ABC(i,:) = balls_to_vols(B2e);
        plot(VF(1),ABC(i,2),'or',VF(1),ABC(i,3),'om')
        i = i+1;
        %disp(i)
        %disp(test)
        %disp(B1)
        %disp(B1e)
        %disp(balls_to_vols(B1e))
        %disp(B2)
        %disp(B2e)
        %disp(balls_to_vols(B2e))
        %disp('__________')
        test = test+1;
    end

end

xlabel('overall vol frac solvent')
ylabel('vol frac P3HT in each phase')
hold off

A = ABC(:,1);
B = ABC(:,2);
C = ABC(:,3);

end

function [B1e,B2e] = equilibrate_bags(B1,B2)

%equilibrate bags takes two 1x3 vectors of integer #s of balls in two bags,
%and moves balls between those bags until the overall energy of the system
%is minimized according to a ?G function

%find the direction of steepest energetic descent, allowing only direction
%which are possible.
%if it is a minima, return the bags
%if not, move balls to a lower energy state and repeat

start = Gibbs_calc(B1,B2,0,0); %make the 1x7 vector of both bags and their starting energy

%create a 7x7 starting matrix to feed to the iterator
initial_grad = [start;zeros(6,7)+1E6]; %add 1,000,000 to the empties so that they don't get mistaken for an energy minima
Eq_bags = iterative_equilibrator(initial_grad); %this is a 1x7 2-bag compositions and energy vector

B1e = Eq_bags(1:3);
B2e = Eq_bags(4:6);

end

function out = iterative_equilibrator(GRAD)

% iterative equilibrator takes a 7x7 GRAD matrix, described below, takes
% the starting composition (first row), generates the six different
% possibilities if a ball were transferred from one bag to the other,
% calculates their energies, and then passes a new 7x7 GRAD matrix to
% itself with the lowest energy compositions as the new start row. If the
% start row is the lowest energy composition, it is returned as 'out' (a
% 1x7 vector)

% [<balls of A in bag 1>... B1B B1C B2A B2B B2C <?G of this composition>]
% first row is the original
% leave all zeros if a direction is impossible

index = 0;

while index ~= 1
    
    for i = 1:3
        if GRAD(1,i)>1 %if there are any balls of i in B1,
            GRAD(2*i,:) = Gibbs_calc(GRAD(1,1:3),GRAD(1,4:6),i,1); %find the energy of moving an i from bag 1
        end
        if GRAD(1,i+3)>1 %if there are any i in B2,
            GRAD(2*i+1,:) = Gibbs_calc(GRAD(1,1:3),GRAD(1,4:6),i,2); %find energy of moving i from bag 2
        end
    end
    
    [lowest_E, index] = min(GRAD(:,7)); %this will output the lowest energy value and its index
    
    if index == 1
        out = GRAD(1,:);
    else
        GRAD = [GRAD(index,:);zeros(6,7)+1E6]; % remake the matrix with the lowest energy as the new starting point
    end
    
end

end

function out = Gibbs_calc(B1,B2,comp,from)

%% Gibbs Calc

% Gibbs calc takes two 1x3 vectors that are the bags, changes the ball
% counts to reflect which *component (scalar, 1 2 or 3) is being moved *from what bag (scalar, 1 or 2),
% calculates the energy of the new conformation, and returns this as a 1x7
% vector:
% [B1A B1B B1C B2A B2B B2C ?G]
% pass a '0' to 'from' to calculate energy without switching balls around

B1new = B1; %start with the new bags the same as given
B2new = B2;

%the following section moves a ball according to which component it is and
%which bag it's taken from

if from == 1
    B1new(comp) = B1(comp) - 1;
    B2new(comp) = B2(comp) + 1;
elseif from == 2
    B1new(comp) = B1(comp) + 1;
    B2new(comp) = B2(comp) - 1;
end

%now we convert total ball counts to mol fractions, then volume fractions
VF1 = balls_to_vols(B1new);
VF2 = balls_to_vols(B2new);

%Energy1 = B1new(1)*log(VF1(1)) + B1new(2)*log(VF1(2)) + B1new(3)*log(VF1(3)) + g12()*B1new(1)*VF1(2) + g13()*B1new(1)*VF1(3) + g23()*B1new(2)*VF1(3);
%Energy2 = B2new(1)*log(VF2(1)) + B2new(2)*log(VF2(2)) + B2new(3)*log(VF2(3)) + g12()*B2new(1)*VF2(2) + g13()*B2new(1)*VF2(3) + g23()*B2new(2)*VF2(3);
%Energy = Energy1+Energy2;
%disp(Energy)

% actual formula: ?G/RT = n1*log(vf1) + n2*log(vf2) + n3*log(vf3) +
% g12(u2)*n1*vf2 + g13(vf3)*n1*vf3 + g23(vf3)*n2vf3
% in this case we are using the ball counts as number of moles. Dividing
% by total number of moles would yield molar gibbs but we're just going to add
% it all up because I'm lazy

R = 8.314;
T = 295;
MV = generate_molar_volumes();
Energy1 = R*T*(B1new(1)*log(VF1(1)) + B1new(2)*log(VF1(2)) + B1new(3)*log(VF1(3)) + (g12()*VF1(1)*VF1(2) + g13()*VF1(1)*VF1(3) + g23()*VF1(2)*VF1(3))*(MV(1)*B1new(1)+MV(2)*B1new(2)+MV(3)*B1new(3)));
Energy2 = R*T*(B2new(1)*log(VF2(1)) + B2new(2)*log(VF2(2)) + B2new(3)*log(VF2(3)) + (g12()*VF2(1)*VF2(2) + g13()*VF2(1)*VF2(3) + g23()*VF2(2)*VF2(3))*(MV(1)*B2new(1)+MV(2)*B2new(2)+MV(3)*B2new(3)));
Energy = Energy1 + Energy2;

% new formula, Prausnitz:
% G = RT*(n1*log(vf1) + n2*log(vf2) + n3*log(vf3) + (X12vf1vf2 + X13vf1vf3
% + X23vf2vf3)*(m1n1 + m2n2 + m3n3))
% m's are molar volume ratios, where m_solvent = 1 (see
% generate_molar_volumes)

out = [B1new, B2new, Energy];

end

function VF = balls_to_vols(B)

molar_vols = generate_molar_volumes();

denom = 0;
for i = 1:3
    denom = denom + B(i)*molar_vols(i);
end

for i = 1:2
    VF(i) = B(i)*molar_vols(i)/denom;
end

VF(3) = 1 - VF(1) - VF(2);

end

function out = g12()
out = Xij('CHCl3','P3HT',295);
%out = 0.99;
end

function out = g13()
out = Xij('CHCl3','PS',295);
%out = 0.39;
end

function out = g23()
out = 0.48;
end

function [B1,B2] = fill_bags(volfracs,balls)

B1 = zeros(1,3); %each bag will have a certain number of balls

comps = vol_to_mol(volfracs); %the starting compositions are volume fractions and must be converted to mol fractions here to fill the bags
B1(1) = ceil(balls*comps(1)); %multiply mol fractions by total balls in each bag
B1(2) = floor(balls*comps(2));
B1(3) = balls - B1(1) - B1(2);

B2 = B1; %both bags start with the same number of balls

end

function comps = vol_to_mol(vf)
%takes a 1x3 vector and outputs another 1x3 vector

molar_vols = generate_molar_volumes(); %molar volumes need not be absolute, just relative [non-solvent solvent polymer]

denom = 0; %first calculate the denominator of the formula
for i = 1:3
    denom = denom + vf(i)/molar_vols(i);
end

for i = 1:2
    comps(i) = (vf(i)/molar_vols(i))/denom;
end

comps(3) = 1 - comps(1) - comps(2);

end

function out = generate_molar_volumes()
% 3 component vector of relative molar volume values
out = [1 119 1635];
end

function out = populate_triangle(steps)

%% Populate Triangle
%this takes a number of steps and outputs an nx3 matrix of compositions
%which are volume fractions along a line of solvent evaporation

out = zeros(steps,3);

pol_vol_fracs = [.0087 .00952]; % vol. frac. P3HT / PS based on 10 mg/mL each
solv_frac = 1 - pol_vol_fracs(1) - pol_vol_fracs(2);

solv_line = linspace(0,solv_frac,steps)';
out(:,1) = solv_line;

for i = 1:steps
    out(i,3) = (1-out(i,1))/(1+pol_vol_fracs(1)/pol_vol_fracs(2));
    out(i,2) = 1 - out(i,1) - out(i,3);
end

end 