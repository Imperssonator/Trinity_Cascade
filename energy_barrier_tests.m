function out = energy_barrier_tests(components,DOP,start_comp,final_comp,steps,Temps)

%% Phase Walk
% This function tracks the thermodynamic driving force and equilibrium
% compositions of phase separation for a process where one component (the solvent) is
% evaporating.
% Components: a cell array containing a string of the name of the compounds
% in the system
% DOPs: the DOP's of each component. 1 = solvent
% Start_comp: A vector containing volume fractions of the corresponding
% starting components
% Final_comp: A vector containing the expected final composition of the
% system, again volume fractions
% Steps: How finely do you want to discretize the process in time?
% balls: how finely do you want to discretize the phases?
% Temps: A vector of length = steps with the temperatures in Kelvin that
% will be experienced by the system at each step. If this vector is not all
% the same constant the resulting phase diagram will not be valid at a
% single temperature

%% Data structure
% [temperature, starting solution vol. fracs, phase 1, phase 2 if applicable]
% [T SVF1 SVF2 SVF3, P1VF1 P1VF2 P1VF3, P2VF1 P2VF2 P2VF3]
% each row is for a different starting solution
% then there's a second matrix that has all of that information as
% discretized balls, i.e. number of moles

global DOPs
DOPs = DOP;

Walk = ntests(start_comp,final_comp,steps,Temps);

Walk = run_tests(Walk);

plot_process(Walk);

out = Walk;

end

function out = plot_process(Walk)

P1 = []; P2 = []; S = []; count = 1;

figure
hold on
for i = 1:length(Walk)
    if Walk(i,5:7,1) == [0 0 0]
        plot(Walk(i,2,1),Walk(i,3,1),'ob',Walk(i,2,1),Walk(i,4,1),'oc')
    else
        plot(Walk(i,2,1),Walk(i,6,1),'or',Walk(i,2,1),Walk(i,9,1),'or',Walk(i,2,1),Walk(i,7,1),'om',Walk(i,2,1),Walk(i,10,1),'om')
        S(count) = Walk(i,5,1);
        P1(count) = Walk(i,6,1);
        P2(count) = Walk(i,7,1);
        count = count+1;
        S(count) = Walk(i,8,1);
        P1(count) = Walk(i,9,1);
        P2(count) = Walk(i,10,1);
    end
end
hold off

figure
ternplot(P1,S,P2,'or')
ternlabel('P3HT','PS','CHCl3')

end

function out = run_tests(Walk)

%% Run Tests
% This will now take our overall compositions and equilibrate, starting
% with the dilute solution. Once phase separation is predicted, it will try
% to find the next point on the binodal from the previous.

% first we equilibrate the first test solution (the most dilute)

for i = 1:length(Walk)
    if i==1 || not(any(Walk(i-1,5:7,1)))
        Walk = single_eqm(Walk,i);
    else
        Walk = guided_eqm(Walk,i);
        disp('ran guided')
        disp(i)
    end
end

out = Walk;
end

function out = single_eqm(Walk,i)

%% Single Equilibrium
% The equilibrator that runs naive, i.e. assumes that its solution is
% already at equilibrium or very close, because the previous solution was
% at equilibrium

T = Walk(i,1,1);
GRAD = [Gibbs_calc(T,Walk(i,2:4,2),Walk(i,2:4,2),0,0); ...
    zeros(6,7)];
% The first GRAD takes the balls in the starting solution at step i and
% doubles them to create two bags (phases) to equilibrate
% Add 1E8 so that the other options aren't mistaken for energy minima

GRAD = iterative_equilibrator(GRAD,T);

if GRAD(1:3) == Walk(i,2:4,2) % if equilibration changed nothing,
    out = Walk; % leave zeros!!! Zeros in columns 5:10 of Walk mean the overall phase is equilibrated
else
    Walk(i,5:10,2) = GRAD(1,1:6);
    Walk(i,5:7,1) = balls_to_vols(Walk(i,5:7,2));
    Walk(i,8:10,1) = balls_to_vols(Walk(i,8:10,2));
    out = Walk;
end
end

function out = guided_eqm(Walk,i)

%% Guided Equilibrium
% This equilibrator runs if the previous test solution phase separated. It
% adds or removes balls equally from both phases until the overall
% composition matches the new test solution, then equilibrates from there,
% which should already be very close. Should.

T = Walk(i,1,1);
Diff = 2*Walk(i,2:4,2) - (Walk(i-1,5:7,2)+Walk(i-1,8:10,2));
P1Diff = ceil(Diff./2);
P2Diff = floor(Diff./2);

GRAD = [Gibbs_calc(T,Walk(i-1,5:7,2)+P1Diff,Walk(i-1,8:10,2)+P2Diff,0,0); ...
    zeros(6,7)];
% Starting GRAD adds the calculated differences between the previous
% binodal points and the new overall composition (most likely just removing
% solvent from both phases) and finds an energy. Then we run iterative eq
% on that composition and hopefully it takes way less iterations.

GRAD = iterative_equilibrator(GRAD,T);

if max(abs(GRAD(1,1:3)-GRAD(1,4:6)))<=1
    out = Walk; % leave zeros because that point is back at eqm
else
    Walk(i,5:10,2) = GRAD(1,1:6);
    Walk(i,5:7,1) = balls_to_vols(Walk(i,5:7,2));
    Walk(i,8:10,1) = balls_to_vols(Walk(i,8:10,2));
    out = Walk;
end

end


function out = iterative_equilibrator(GRAD,T)

%% Iterative Equilibrator
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

if find(GRAD(1,1:6)==0)
    out = GRAD;
    return
end

while index ~= 1
    
    for i = 1:3
        if GRAD(1,i)>1 %if there are any balls of i in B1,
            GRAD(2*i,:) = Gibbs_calc(T,GRAD(1,1:3),GRAD(1,4:6),i,1); %find the energy of moving an i from bag 1
            if GRAD(2*i,7)<GRAD(1,7)
                disp('moving...')
                disp(GRAD)
            end
        end
        if GRAD(1,i+3)>1 %if there are any i in B2,
            GRAD(2*i+1,:) = Gibbs_calc(T,GRAD(1,1:3),GRAD(1,4:6),i,2); %find energy of moving i from bag 2
            if GRAD(2*i+1,7)<GRAD(1,7)
                disp('moving...')
                disp(GRAD)
            end
        end
    end
    
    %disp(GRAD)
    [lowest_E, index] = min(GRAD(:,7)); %this will output the lowest energy value and its index
    
    if index == 1
        out = GRAD;
    else
        GRAD = [GRAD(index,:);zeros(6,7)+1E6]; % remake the matrix with the lowest energy as the new starting point
    end
    
end

end

function out = Gibbs_calc(T,B1,B2,comp,from)

%% Gibbs Calc

% Gibbs calc takes two 1x3 vectors that are the bags, changes the ball
% counts to reflect which *component (scalar, 1 2 or 3) is being moved *from what bag (scalar, 1 or 2),
% calculates the energy of the new conformation, and returns this as a 1x7
% vector:
% [B1A B1B B1C B2A B2B B2C ?G]
% pass a '0' to 'from' to calculate energy without switching balls around

global DOPs
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

R = 8.314;
k = 1.38E-23;
T = 273;

Energy1 = B1new(1)*log(VF1(1)) + B1new(2)*log(VF1(2)) + B1new(3)*log(VF1(3)) + g12()*B1new(1)*VF1(2) + g13()*B1new(1)*VF1(3) + g23()*B1new(2)*VF1(3);
Energy2 = B2new(1)*log(VF2(1)) + B2new(2)*log(VF2(2)) + B2new(3)*log(VF2(3)) + g12()*B2new(1)*VF2(2) + g13()*B2new(1)*VF2(3) + g23()*B2new(2)*VF2(3);
Energy = R*T*(Energy1+Energy2);
%disp(Energy)

% actual formula: ?G/RT = n1*log(vf1) + n2*log(vf2) + n3*log(vf3) +
% g12(u2)*n1*vf2 + g13(vf3)*n1*vf3 + g23(vf3)*n2vf3
% in this case we are using the ball counts as number of moles. Dividing
% by total number of moles would yield molar gibbs but we're just going to add
% it all up because I'm lazy


%Energy1 = (B1new(1)*log(VF1(1)) + B1new(2)*log(VF1(2)) + B1new(3)*log(VF1(3)) + (g12()*VF1(1)*VF1(2) + g13()*VF1(1)*VF1(3) + g23()*VF1(2)*VF1(3))*(DOPs(1)*B1new(1)+DOPs(2)*B1new(2)+DOPs(3)*B1new(3)));
%Energy2 = (B2new(1)*log(VF2(1)) + B2new(2)*log(VF2(2)) + B2new(3)*log(VF2(3)) + (g12()*VF2(1)*VF2(2) + g13()*VF2(1)*VF2(3) + g23()*VF2(2)*VF2(3))*(DOPs(1)*B2new(1)+DOPs(2)*B2new(2)+DOPs(3)*B2new(3)));
%Energy = R*T*(Energy1 + Energy2);
%disp(Energy)

% new formula, Prausnitz:
% G = RT*(n1*log(vf1) + n2*log(vf2) + n3*log(vf3) + (X12vf1vf2 + X13vf1vf3
% + X23vf2vf3)*(m1n1 + m2n2 + m3n3))
% m's are molar volume ratios, where m_solvent = 1 (see
% generate_molar_volumes)

out = [B1new, B2new, Energy];

end

function out = g12()
%out = Xij('CHCl3','P3HT',295);
out = 0.99; %CHCl3/P3HT
%out = 0.4;
end

function out = g13()
%out = Xij('CHCl3','PS',295);
out = 0.39; % CHCl3/PS
%out = 0.4;
end

function out = g23()
out = 0.48; %PS/P3HT
%out = 0.004;
end

function out = ntests(start,final,steps,Temps)
%% ntests
% this takes starting and final compositions and creates a discretized
% linear path between the two in terms of volume fractions

out = zeros(steps,10,2);
out(:,1,1) = Temps;
out(1,2:4,1) = start;
out(steps,2:4,1) = final;

path_vector = final - start; %find a linear path in component vol. frac. space between start and finish

for i = 2:steps-1
    out(i,2:4,1) = start+path_vector*(i-1)/(steps-1);
end

out = vols_to_balls(out);

end 

function out = vols_to_balls(Walk)
%% Vols to Balls
% input: the Walk matrix with the vol. frac. path already initialized
% output: the Walk matrix with ball counts calculated from vol. frac's

global DOPs
MFstart = vols_to_mols(Walk(1,2:4,1));
MFend = vols_to_mols(Walk(end,2:4,1));

% The basis is the index of the component that is both a polymer and has
% the lowest mole fraction at the end of the process
polymers = find(DOPs>1);
[foo basis] = min(MFend(polymers));
if length(polymers)==2
    basis = basis+1;
else
    basis = basis+2;
end

Walk(end,2:4,2) = mols_to_balls(MFend,basis);

for i = 1:length(Walk)-1
    Walk(i,2:4,2) = mols_to_balls(vols_to_mols(Walk(i,2:4,1)),basis);
end

out = Walk;

end

function out = mols_to_balls(MF,basis)

%% Mols to Balls
% take mol fraction vector and the component that is the basis and spit out
% a vector of ball counts with 100 for the basis component
base_count = 10;

for i = 1:3
    if i == basis
        out(i) = base_count;
    else
        out(i) = round(MF(i)/MF(basis)*base_count);
    end
end
end

function mol_fracs = vols_to_mols(vf)

%% Vols to mols
% volume fractions converted to mol fractions

global DOPs
denom = 0; %first calculate the denominator of the formula
mol_fracs = zeros(1,3);
for i = 1:3
    denom = denom + vf(i)/DOPs(i);
end

for i = 1:2
    mol_fracs(i) = (vf(i)/DOPs(i))/denom;
end

mol_fracs(3) = 1 - mol_fracs(1) - mol_fracs(2);

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