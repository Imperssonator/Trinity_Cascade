function out = E_barrier_tests()



Temp=301;
species = {'P3HT'; 'CHCl3'}; %chemical species in system
DP = [236; 1]; %relative molar volumes of each species... indicating P3HT has a DoP of 236
VFO = [0.01; 0.99]; %We will change this on every iteration of the for loop coming up
nuc = 1E-12;
x0 = [1-nuc nuc;1-nuc nuc];
PFeq = [.2535 .7465; 0.9745 0.0255];
vec = PFeq-x0
X = Chi(species,Temp);

pts = linspace(0,0.000001,1000);

figure
hold on
for i = pts
    plot(i,GIBBS2bin(x0+vec*i,VFO,DP,X),'-b')
end

% figure
% hold on
% h=1000;
% pts = linspace(0,1,h);
% disp(X)
% En = [];
% for i = pts
%     VF = [i; 1-i];
%     E = SINGLE_GIBBS(VF,DP,X);
%     plot(i,E(1),'-b',i,E(2),'-b',i,E(3),'-r');
%     En = [En,E];
% end