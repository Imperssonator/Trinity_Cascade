%disp('2 polymer 1 solvent')
%W = phase_walk([0 0 0],[1 1000 1000],[0.99 0.005 0.005],[0 0.5 0.5],100,zeros(100,1)+273);
disp('CHCl3 P3HT PS')
W = phase_walk([0 0 0],[1 119 1635],[0.99 0.005 0.005],[0 0.5 0.5],100,zeros(100,1)+292);