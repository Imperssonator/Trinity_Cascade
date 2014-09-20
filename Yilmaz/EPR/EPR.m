function [VFeq,PFeq,Eeq,stab,VFspin] = EPR(system)
%% Enthalpic Perturbation/Relaxation
% Finds global minima in the Gibbs function by strengthening enthalpic
% contributions until metastable states become unstable, then uses the
% separated phase as the new initial guess as the enthalpic terms relax
% back to their initial values.
%
% DP is 3X1! Important!
% components are always rows
% phases are always columns
% 
% system is a structure. system.VFO is the overall volume fraction, a 3x1
% system.DP is the 3x1 of degrees of polymerization (relative molar
% volumes)
% system.species is a 3x1 cell array of strings with the species names

x0 = system.x0;
DP = system.DP;
VFO = system.VFO;
species = system.species;
T = system.Temp;

if length(species) == 3
    [VFeq,PFeq,Eeq,stab] = binodal_epr(DP,VFO,x0,species);
elseif length(species) == 2
    [VFeq,PFeq,Eeq,stab] = binary_epr(DP,VFO,x0,species,T);
end

% if stab == 1 || stab == 0
%     VFspin = spinodal(DP,VFO,PFeq,species);
% else
%     VFspin = zeros(3,2);
% end

% VFcrit = criticalpt(DP);
% disp('critical point:')
% disp(VFcrit)

end

function [VFeq,PFeq,Eeq,stab] = binodal_epr(DP,VFO,x0,species)
%% Binodal perturbation relaxation
% Start by determining the energy of a single phase
% Then, determine stability by running fmincon.
% If unstable, binodal should be determined instantly... heh heh heh...
% If stable, increase the pert and try again, until you reach the pert
% limit, then it must be stable.
% If you reach an unstable phase before reaching the pert limit, find its
% eqm phase fractions, then use those as starting guesses as you relax back
% to 0 perturbation. Hopefully it will be stable and get you to the global
% energy minimum.

[PFi,Ei,stabi] = SUPERGRADSEARCH(x0,VFO,DP,0,.001,1E-6,species);
if stabi == 0
    %disp(VFO)
    disp('was unstable')
    stab = 0;
    PFeq = PFi;
    VFeq = PF2VF(PFi,VFO);
    Eeq = Ei;
else
    disp('running EPR')
    [PFeq,Eeq,stab] = epr_iter(DP,VFO,x0,Ei,species);
    VFeq = PF2VF(PFeq,VFO);
end

end

function [VFeq,PFeq,Eeq,stab] = binary_epr(DP,VFO,x0,species,T)
%% Binodal perturbation relaxation
% Start by determining the energy of a single phase
% Then, determine stability by running fmincon.
% If unstable, binodal should be determined instantly... heh heh heh...
% If stable, increase the pert and try again, until you reach the pert
% limit, then it must be stable.
% If you reach an unstable phase before reaching the pert limit, find its
% eqm phase fractions, then use those as starting guesses as you relax back
% to 0 perturbation. Hopefully it will be stable and get you to the global
% energy minimum.

[PFi,Ei,stabi] = GRADSEARCH2(zeros(2,2)+0.5,VFO,DP,0,.001,1E-10,species,T); %do the first grad search with 0.5's as x0 otherwise it won't predict metastability ever
if stabi == 0
    %disp(VFO)
    disp('was unstable')
    stab = 0;
    [PFeq,Eeq,stab] = GRADSEARCH2(x0,VFO,DP,0,.001,1E-10,species,T); %do the first grad search with 0.5's as x0 otherwise it won't predict metastability ever
    VFeq = PF2VF(PFeq,VFO);
%     PFeq = PFi;
%     VFeq = PF2VF(PFi,VFO);
%     Eeq = Ei;
else
    disp('running EPR')
    [PFeq,Eeq,stab] = epr_iter_bin(DP,VFO,x0,Ei,species,T); %but if you must do an EPR, use the previous point's PFeq as x0 (fed in the system struct)
    VFeq = PF2VF(PFeq,VFO);
end

end

function [PFeq,Eeq,stab] = epr_iter(DP,VFO,x0,Ei,species)
%% EPR Iterative Search
% Implements enthalpic perturbations, uses super grad search to find minima
% until phase separation is predicted, then relaxes back down, re-locating
% minima on the way

pert = 0; maxpert = 3; pertstep = 0.001; stab = 1; difftol = 1E-6; initstep = .001;
conTol = 1E-6;

disp('perturbing...')
while stab > 0 && pert<=maxpert
    pert = pert+pertstep;
    [PFiter,Eiter,stab] = SUPERGRADSEARCH(x0,VFO,DP,pert,initstep,conTol,species);
end

if pert>=maxpert
    PFeq = zeros(3,2)+0.5;
    Eeq = Ei;
    stab = 2;
    disp('was absolutely stable')
    return
end

initstep = 1E-4; % could change these values for the trip back down
conTol = 1E-7;

disp('relaxing...')
while pert>0
    pert = pert-pertstep;
    if pert<=0.01
        %disp('almost done')
        conTol = 1E-7;
    end
    x0 = PFiter;
    [PFiter,Eiter,stab] = SUPERGRADSEARCH(x0,VFO,DP,pert,initstep,conTol,species);
end

if PFDIFF(PFiter)>difftol && Eiter<Ei
    stab = 1;
    disp('was metastable!')
    PFeq = PFiter;
    Eeq = Eiter;
else
    stab = 2;
    disp('was absolutely stable')
    VFeq = [VFO,VFO];
    PFeq = zeros(3,2)+0.5;
    Eeq = Ei;
end
end

function [PFeq,Eeq,stab] = epr_iter_bin(DP,VFO,x0,Ei,species,T)
%% EPR Iterative Search
% Implements enthalpic perturbations, uses super grad search to find minima
% until phase separation is predicted, then relaxes back down, re-locating
% minima on the way

pert = 0; maxpert = 3; pertstep = 0.001; stab = 1; difftol = 1E-6; initstep = .001;
conTol = 1E-10;
disp('perturbing...')
while stab > 0 && pert<=maxpert
    pert = pert+pertstep;
    [PFiter,Eiter,stab] = GRADSEARCH2(x0,VFO,DP,pert,initstep,conTol,species,T);
end

if pert>=maxpert
    PFeq = zeros(2,2)+0.5;
    Eeq = Ei;
    stab = 2;
    disp('was absolutely stable')
    return
end

initstep = 1E-6; % could change these values for the trip back down
conTol = 1E-10;

disp('relaxing...')
while pert>0
    pert = pert-pertstep;
    if pert<=0.01
        %disp('almost done')
        conTol = 1E-10;
    end
    x0 = PFiter;
    [PFiter,Eiter,stab] = GRADSEARCH2(x0,VFO,DP,pert,initstep,conTol,species,T);
end

if abs(PFiter(1,1)-PFiter(2,1))^2>difftol && Eiter<Ei % if phase I actually has different fractions of each component (i.e. isn't the same composition as Phase II but of different size) and the system actually has lower energy...  
    stab = 1; % then it was metastable.
    disp('was metastable!')
    PFeq = PFiter;
    Eeq = Eiter;
else
    stab = 2;
    disp('was absolutely stable')
    VFeq = [VFO,VFO];
    PFeq = zeros(2,2)+0.5;
    Eeq = Ei;
end
end

function out = spinodal(DP,VFO,x0,species)
%% Spinodal
% If VFO was unstable or metastable, use the binodal PFeq as an initial
% guess for the spinodal curve. Try PF's along a vector that points from
% the binodal inward towards VFO, because the solution must be along that
% line



% options = optimoptions('fmincon','MaxFunEvals',10000,'Algorithm','interior-point','TolCon',1E-8);
% S_obj = @(x) SPINOBJ(x,VFO,DP);
% lb = zeros(3,1)+1E-30;
% ub = zeros(3,1)+1-1E-30;
% PROBLEM = createOptimProblem('fmincon','objective',S_obj,'x0',x0,'lb',lb,'ub',ub,'options',options);
% GS = GlobalSearch;
% [PFspin,ShouldBeZero] = run(GS,PROBLEM);
% PFspin = [PFspin,1-PFspin];
% 
% disp(ShouldBeZero)
% disp('should be zero')
% disp(PFspin)
% VFspin = PF2VF(PFspin,VFO);
% out = VFspin;

end


function out = criticalpt(DP)
%% Critical Point
% This will probably never do anything

%x0 = [0.8; 0.1; 0.1]; gets it for Yilmaz
x0 = [0.25; 0.25; 0.5];
options = optimoptions('fmincon','MaxFunEvals',10000,'Algorithm','interior-point','TolCon',1E-8);
C_obj = @(x) CRITOBJ(x,DP);
lb = zeros(3,1)+1E-30;
ub = zeros(3,1)+1-1E-30;
PROBLEM = createOptimProblem('fmincon','objective',C_obj,'x0',x0,'lb',lb,'ub',ub,'options',options);
GS = GlobalSearch;
[VFcrit,fval] = run(GS,PROBLEM);
disp('should be zero:')
disp(fval)
out = VFcrit;

end