function [dE,transfer] = LOCGRAD(VFO,PF,DP,dir,step)
%% LOCGRAD
% VFO is 1x3 [overall vol. fracs of each component]
% PF is 3x2 phase fractions
% DP is "degrees of polymerization".. 1x3
%
% dir is a 3x1 vector indicating the direction of movement in phase space.
% [1;0;0] means moving component 1 from phase 1 into phase 2,
% [1;-1;0.5] means move equal volumes comp. 1 into phase 2 and comp. 2 into
% phase 1, and move half that much comp 3 into phase 2.
%
% step indicates how large of a step you'd like to take in that direction.
% step will be altered if it would take PF into negative numbers.
%
% out is the energy difference cause by that move

norm_dir = dir/norm(dir); %normalize the input vector
step_vec = step*norm_dir; %create the step vector


N_current = PF2N(PF,VFO,DP);
PF_new = 

end