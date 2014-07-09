function out = phase_walk(components,DOP,start_comp,final_comp,steps,Temps)

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
% the second layer of Walk works in terms of moles:
% [Energy, [total moles of each compound in system], [molar fraction of each
% compound in PHASE 1]
% so the second layer only uses out to its 7th index, where the first layer
% uses out to 10

global DOPs
DOPs = DOP;

Walk = ntests(start_comp,final_comp,steps,Temps);

out=TL_Energy(Walk,75);
return

Walk = run_tests(Walk);

plot_process(Walk);

out = Walk;

end

function out = plot_process(Walk)

P1 = []; P2 = []; S = []; count = 1;

figure
hold on
for i = 1:length(Walk)
    if not(any(Walk(i,5:7,1)))
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
    Walk = single_eqm(Walk,i);
end

% for i = 1:length(Walk)
%     if i==1 || not(any(Walk(i-1,5:7,1)))
%         Walk = single_eqm(Walk,i);
%     else
%         Walk = guided_eqm(Walk,i);
%         disp('ran guided')
%         disp(i)
%     end
% end

out = Walk;
end

function out = TL_Energy(Walk,i)
%% TL Energy
% Trying to visualize the energy as a function of degree of phase
% separation

global DOPs
T = Walk(i,1,1);
start_mols = Mtot(Walk,i);
if all(start_mols)
    figure
    hold on
    for i = 0:0.01:0.1
        X = [0.5 0.5+i*0.5 0.5-i*0.5];
        Energy=Gibbs_calc(X,start_mols,DOPs,T);
        plot(i,Energy,'or')
    end 
end
out = 1;
end

function out = single_eqm(Walk,i)

%% Single Equilibrium
% The equilibrator that runs naive, i.e. assumes that its solution is
% already at equilibrium or very close, because the previous solution was
% at equilibrium

global DOPs
T = Walk(i,1,1);
start_mols = Mtot(Walk,i);
if all(start_mols)
    f = @(x)Gibbs_calc(x,start_mols,DOPs,T);
    %options = saoptimset('TolFun',1e-8);
    [eqm_mfs,Energy] = fmincon(f,[0.5 0.5 0.5],[],[],[],[],[0 0 0],[1 1 1],[],optimoptions('fmincon','Algorithm','interior-point'));
    %[eqm_mfs,Energy] = simulannealbnd(f,[0.5 0.5 0.5],[0 0 0],[1 1 1]);
    Walk(i,1,2) = Energy;
    %disp(Energy)
    if eqm_mfs ~= [0.5 0.5 0.5] % if the eqm mfs change from start, store the new data in Walk
        Walk(i,5:7,2) = eqm_mfs;
        Walk(i,5:10,1) = mols_to_vols(Walk(i,2:7,2));
    end
else
    disp(['component ' num2str(find(start_mols==0)) 'is empty'])
end

out = Walk;

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

function out = ntests(start,final,steps,Temps)
%% ntests
% this takes starting and final compositions and creates a discretized
% linear path between the two in terms of volume fractions
% it also generates the Walk matrix

out = zeros(steps,10,2);
out(:,1,1) = Temps;
out(1,2:4,1) = start;
out(steps,2:4,1) = final;

path_vector = final - start; %find a linear path in component vol. frac. space between start and finish

for i = 2:steps-1
    out(i,2:4,1) = start+path_vector*(i-1)/(steps-1);
end

out = make_mols(out);

end 

function out = make_mols(Walk)
%% Make Mols
% input: the Walk matrix with the vol. frac. path already initialized
% output: the Walk matrix with mole counts calculated from vol. frac's,
% specifically [0, total moles of each comp., zeros....]

Mtotal = 1; % hard-coded total number of moles in system... shouldn't really matter

for i = 1:length(Walk)
    MFi = vols_to_mfs(Walk(i,2:4,1)); % convert vol fracs to mol fracs, then to total mols with a system total of 100
    Walk(i,2:3,2) = MFi(1:2).*100;
    Walk(i,4,2) = 100 - sum(Walk(i,2:3,2));
end

out = Walk;

end

function mol_fracs = vols_to_mfs(vf)

%% Vols to mfs
% volume fractions converted to mol fractions, 1x3 vectors only

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

function out = mols_to_vols(M)

%% mols to vols
% takes a 1x6 [total moles in system, molar fraction of each comp. in phase
% 1] and outputs another 1x6, [VFs phase 1, VFs phase 2]

M1 = M(1:3).*M(4:6);
M2 = M(1:3) - M1;

out(1:3) = M2V(M1);
out(4:6) = M2V(M2);

end

function out = M2V(M)
% just an inner function of mols to vols that needs to be run twice,
% converts 1x3 vector of moles to 1x3 vector of VFs

global DOPs
denom = 0;
for i = 1:3
    denom = denom+M(i)*DOPs(i);
end
for i = 1:2
    out(i) = M(i)*DOPs(i)/denom;
end
out(3) = 1 - out(1) - out(2);
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
%out = 0.99; %CHCl3/P3HT
out = 0.4;
end

function out = g13()
%out = Xij('CHCl3','PS',295);
%out = 0.39; % CHCl3/PS
out = 0.44;
end

function out = g23()
%out = 0.48; %PS/P3HT
out = 0.004;
end

function out = Mtot(Walk,i)
% get total number of moles of each component at point i
out = Walk(i,2:4,2);

end

function out = Vtot(Walk,i)
% get overall volume fractions of each comp. at point i
out = Walk(i,2:4,1);
end

function out = VPi(Walk,test,i)
% get vol. fracs. of each comp in phase i
if i == 1
    out = Walk(test,5:7,1);
else
    out = Walk(test,8:10,1);
end
end
