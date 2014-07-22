function [VFeq,VFspin] = Yilmaz(DP,species)
%% Yilmaz
% Implements Yilmaz's equations for finding binodal and spinodal lines
% DP is 3X1! Important!
% components are always rows
% phases are always columns

VFO = [0.3 0.3 0.4]'; % Yilmaz Fig 1a or Prausnitz 1b... unstable in both
x0 = [0.5 0.5 0.5]'; % good init guess for Praus 1b

[VFeq,PFeq] = binodal(DP,VFO,x0,species);
% VFspin = spinodal(DP,VFO,PFeq(:,1));
% VFcrit = criticalpt(DP);
% disp('critical point:')
% disp(VFcrit)

% SPIN = spinodal(DOP,T,species);
% CRIT = critpt(DOP,T,species);

end

function [VFeq,PFeq] = binodal(DP,VFO,x0,species)
%% Binodal

options = optimoptions('fmincon','MaxFunEvals',10000,'Algorithm','interior-point','TolCon',1E-8);
B_obj = @(x) BINOBJ(x,VFO,DP);
lb = zeros(3,1)+1E-30;
ub = zeros(3,1)+1-1E-30;
PROBLEM = createOptimProblem('fmincon','objective',B_obj,'x0',x0,'lb',lb,'ub',ub,'options',options);
GS = GlobalSearch('NumStageOnePoints',200,'NumTrialPoints',1000);
[PFeq,Energy] = run(GS,PROBLEM);
PFeq = [PFeq,1-PFeq];
disp(PFeq)
disp(Energy)
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