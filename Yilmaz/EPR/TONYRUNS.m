function out = TONYRUNS()

%% Build the structure input

system1 = struct();
system1.species = {'CHCl3'; 'hexane'; 'P3HT'}; %chemical species in system
system1.DP = [1; 1; 236]; %relative molar volumes of each species... indicating P3HT has a DoP of 236
system1.VFO = [0; 0; 0;]; %We will change this on every iteration of the for loop coming up
points = 50; % # random points to test per phase diagram
save('input_struct.mat','system1')

%% Pick system volume fractions to test
% The "EPR" function mimics reality in that it takes a test solution
% defined by overall volume fractions of each component, a column in VFOs,
% and tries to phase separate it by creating two initial phases that each have
% half the initial volume of each component. It then tries moving small
% amounts of each component back and forth until it finds an energetic gradient
% it can descend until it has two phases whose combined energy is less than
% the starting energy. So to make a whole phase diagram it needs about
% 50-100 test solutions, which we define randomly here. It's like throwing
% darts at the ternary diagram.

VFOs = zeros(3,points); %this will contain all the random overall volume fractions we're going to throw at the program

for i = 1:points
    r1 = rand; r2 = rand;
    VFOs(1,i) = r1; % random vf1
    VFOs(2,i) = (1-r1)*r2; % use r2 as the fraction of the remaining volume after vf1 has been set
    VFOs(3,i) = 1-VFOs(1,i)-VFOs(2,i); % vf3 is the balance
end

%% Build the phase diagram
% run all the test systems. If they phase separate, store the final eqm.
% volume fractions, if not, don't store.

TP = []; %this is the "out". It will contain 3x1 vectors that represent eqm. volume fractions that fall along the binodal curve. Each "VFO" could produce two vectors in TP

for i = 1:points
    disp('_______________')
    disp(i)
    disp(VFOs(:,i))
    system1.VFO = VFOs(:,i);
    [VFeq,PFeq,Eeq,stab] = EPR(system1);
    disp(VFeq)
    %disp(PFeq)
    if stab == 0
        %plot(x,VFeq(1,1),'or',x,VFeq(1,2),'xr',x,VFeq(2,1),'ok',x,VFeq(2,2),'xk','MarkerSize',12)
        TP = [TP,VFeq];
    elseif stab == 1
        %plot(x,VFeq(1,1),'ob',x,VFeq(1,2),'xb',x,VFeq(2,1),'oc',x,VFeq(2,2),'xc','MarkerSize',12)
        TP = [TP,VFeq];
    else
        %plot(x,VFeq(1,1),'og',x,VFeq(2,1),'oy','MarkerSize',12)
    end
end

out = TP;
save('out.mat','out')
figure
ternplot(TP(1,:),TP(3,:),TP(2,:),'ob','MarkerSize',10) %plot solvent, polymer, non-solvent on ternary diagram
%vertexlabel('hexane','CHCl3','P3HT')

end