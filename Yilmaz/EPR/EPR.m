function [VFeq,PFeq,Eeq,stab] = EPR(VFO,DP,species)
%% Enthalpic Perturbation/Relaxation
% Finds global minima in the Gibbs function by strengthening enthalpic
% contributions until metastable states become unstable, then uses the
% separated phase as the new initial guess as the enthalpic terms relax
% back to their initial values.
%
% DP is 3X1! Important!
% components are always rows
% phases are always columns

x0 = zeros(3,2)+0.5;

[VFeq,PFeq,Eeq,stab] = binodal_epr(DP,VFO,x0,species);
% VFspin = spinodal(DP,VFO,PFeq(:,1));
% VFcrit = criticalpt(DP);
% disp('critical point:')
% disp(VFcrit)

% SPIN = spinodal(DOP,T,species);
% CRIT = critpt(DOP,T,species);

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

[PFi,Ei,stabi] = MYGRADSEARCH(x0,VFO,DP,0,.001,1E-6);
if stabi == 0
    %disp(VFO)
    disp('was unstable')
    stab = 0;
    PFeq = PFi;
    VFeq = PF2VF(PFi,VFO);
    Eeq = Ei;
else
    disp('running EPR')
    [PFeq,Eeq,stab] = epr_iter(DP,VFO,x0,Ei);
    VFeq = PF2VF(PFeq,VFO);
end

end

function [PFeq,Eeq,stab] = epr_iter(DP,VFO,x0,Ei)
%% EPR Iterative Search
% Implements enthalpic perturbations

pert = 0; maxpert = 3; pertstep = 0.001; stab = 1; difftol = 1E-6; initstep = .001;
conTol = 1E-8;

disp('perturbing...')
while stab > 0 && pert<=maxpert
    pert = pert+pertstep;
    [PFiter,Eiter,stab] = MYGRADSEARCH(x0,VFO,DP,pert,initstep,conTol);
end

if pert>=maxpert
    PFeq = zeros(3,2)+0.5;
    Eeq = Ei;
    stab = 2;
    disp('was absolutely stable')
    return
end

initstep = 1E-3;
conTol = 1E-6;

disp('relaxing...')
while pert>0
    pert = pert-pertstep;
    if pert<=0.01
        disp('almost done')
        conTol = 1E-6;
    end
    x0 = PFiter;
    [PFiter,Eiter,stab] = MYGRADSEARCH(x0,VFO,DP,pert,initstep,conTol);
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
    PFeq = x0;
    Eeq = Ei;
end
end

function [VFeq,PFeq,Energy] = binodal_grad(DP,VFO,x0,pert,species)
%% Binodal Gradient Search
% Will only find local minima

% B_obj = @(x) BINOBJ(x,VFO,DP,pert);
% lb = zeros(3,1)+1E-30;
% ub = zeros(3,1)+1-1E-30;
% options = optimoptions('fmincon','MaxFunEvals',10000,'Algorithm','interior-point','TolCon',1E-8);
% PROBLEM = createOptimProblem('fmincon','objective',B_obj,'x0',x0,'lb',lb,'ub',ub,'options',options);
% GS = GlobalSearch('NumStageOnePoints',200,'NumTrialPoints',1000);
% [PFeq,Energy] = run(GS,PROBLEM);
% [PFeq,Energy] = fmincon(PROBLEM);



PFeq = [PFeq,1-PFeq];
% disp(PFeq)
% disp(Energy)
VFeq = PF2VF(PFeq,VFO);

end

function out = spinodal(DP,VFO,x0)
%% Spinodal

options = optimoptions('fmincon','MaxFunEvals',10000,'Algorithm','interior-point','TolCon',1E-8);
S_obj = @(x) SPINOBJ(x,VFO,DP);
lb = zeros(3,1)+1E-30;
ub = zeros(3,1)+1-1E-30;
PROBLEM = createOptimProblem('fmincon','objective',S_obj,'x0',x0,'lb',lb,'ub',ub,'options',options);
GS = GlobalSearch;
[PFspin,ShouldBeZero] = run(GS,PROBLEM);
PFspin = [PFspin,1-PFspin];

disp(ShouldBeZero)
disp('should be zero')
disp(PFspin)
VFspin = PF2VF(PFspin,VFO);
out = VFspin;

end


function out = criticalpt(DP)
%% Critical Point

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