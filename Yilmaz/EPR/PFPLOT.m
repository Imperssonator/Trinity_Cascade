VFO = [.495 0.495 0.01]';
DP = [1000 1000 1]';
G1 =  [0 0.004 0.4;...
       0 0     0.44;...
       0 0     0];
G2 = G1;
pert = 0;
G1 = G1+pert;
G2 = G2+pert;
X = linspace(0.49,0.51,100);
figure
hold on
for i = 1:100
    plot(X(i),GIBBS2([X(i) 1-X(i); 1-X(i) X(i); 0.5 0.5;],VFO,DP,G1,G2),'-b')
end

