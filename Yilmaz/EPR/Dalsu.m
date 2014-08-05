function [comps,DOP] = Dalsu(VP,mg)
%20mg P3HT/2mL CHCl3, diluted down to whatever VP% of 2-MP

DOP = [1 1.64 119]';
comps=zeros(1,3)';
comps(2) = VP/100;
comps(3) = mg/(2/(1-VP/100))/1000/1.1;
comps(1) = 1-comps(2)-comps(3);
end