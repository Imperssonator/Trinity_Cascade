function out = GANGRUNS()

%% Build the structure input

system = struct();
species = {'P3HT'; 'CHCl3'}; %chemical species in system
DP = [236; 1]; %relative molar volumes of each species... indicating P3HT has a DoP of 236
VFO = [0.01; 0.99]; %We will change this on every iteration of the for loop coming up
points = 100; % # number of temperatures to run
% save('gang_input.mat','system1')

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

T_start = 273;
T_end = 373;
Temps = linspace(T_start,T_end,points); %this will contain all the random overall volume fractions we're going to throw at the program

%% Build the phase diagram
% run all the test systems. If they phase separate, store the final eqm.
% volume fractions, if not, don't store.

TP = []; %this is the "out". It will contain 3x1 vectors that represent temperatures and eqm. volume fractions that fall along the binodal curve. Each "VFO" could produce two vectors in TP

for i = 1:points
    disp('_______________')
    disp(i)
    disp(Temps(i))
    system(i).Temp = Temps(i);
    system(i).species = species;
    system(i).DP = DP;
    system(i).VFO = VFO;
    [VFeq,PFeq,Eeq,stab] = EPR(system(i));
    disp(VFeq)
    %disp(PFeq)
    if stab == 0
        %plot(x,VFeq(1,1),'or',x,VFeq(1,2),'xr',x,VFeq(2,1),'ok',x,VFeq(2,2),'xk','MarkerSize',12)
        TP = [TP,[Temps(i),Temps(i); VFeq]];
    elseif stab == 1
        %plot(x,VFeq(1,1),'ob',x,VFeq(1,2),'xb',x,VFeq(2,1),'oc',x,VFeq(2,2),'xc','MarkerSize',12)
        TP = [TP,[Temps(i),Temps(i); VFeq]];
    else
        %plot(x,VFeq(1,1),'og',x,VFeq(2,1),'oy','MarkerSize',12)
    end
end

out = TP;
%save('out.mat','out')
figure
plot(TP(2,:),TP(1,:),'ob')
xlabel('volume fraction P3HT')
ylabel('Temperature (K)')
%vertexlabel('hexane','CHCl3','P3HT')

end