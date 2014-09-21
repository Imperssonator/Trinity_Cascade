function out = GANGRUNS()

%% Build the structure input

system = struct();
species = {'P3HT'; 'CHCl3'}; %chemical species in system
DP = [236; 1]; %relative molar volumes of each species... indicating P3HT has a DoP of 236
VFO = [0.005; 0.995]; %We will change this on every iteration of the for loop coming up
x0 = zeros(2,2)+0.5;
points = 150; % # number of temperatures to run
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

START = 1;
system(START).x0 = x0; % equilibrate the first point from half and half phase split
END = points;

for i = START:END
    disp('_______________')
    disp(i)
    disp(Temps(i))
    disp(VFO)
    disp(system(i).x0)
    system(i).Temp = Temps(i);
    system(i).species = species;
    system(i).DP = DP;
    system(i).VFO = VFO;
    [VFeq,PFeq,Eeq,stab] = EPR(system(i));
    disp(VFeq)
    disp(PFeq)
    if stab == 0
        %plot(x,VFeq(1,1),'or',x,VFeq(1,2),'xr',x,VFeq(2,1),'ok',x,VFeq(2,2),'xk','MarkerSize',12)
        TP = [TP,[Temps(i),Temps(i); VFeq]];
        system(i+1).x0 = PFeq; %use the currently equilibrated phase fractions as initial guess for next temperature
%         system(i+1).x0 = zeros(2)+0.5;
    elseif stab == 1
        %plot(x,VFeq(1,1),'ob',x,VFeq(1,2),'xb',x,VFeq(2,1),'oc',x,VFeq(2,2),'xc','MarkerSize',12)
        TP = [TP,[Temps(i),Temps(i); VFeq]];
        system(i+1).x0 = PFeq; %use the currently equilibrated phase fractions as initial guess for next temperature
%         system(i+1).x0 = zeros(2)+0.5;
    else
        break
        %plot(x,VFeq(1,1),'og',x,VFeq(2,1),'oy','MarkerSize',12)
        system(i+1).x0 = zeros(2,2)+0.5;
    end
end

out = TP;
%save('out.mat','out')
figure
plot(TP(2,:),TP(1,:),'ob','MarkerSize',8)
xlabel('volume fraction P3HT')
ylabel('Temperature (K)')
%vertexlabel('hexane','CHCl3','P3HT')

end