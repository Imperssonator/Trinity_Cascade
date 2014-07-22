function [comps,DOP] = Dalsu(VP,mg)

DOP = [72.7 1 1.64];
comps=zeros(1,3);
comps(3) = VP/100;
comps(1) = mg/(2/(1-VP/100))/1000/1.1;
comps(2) = 1-comps(1)-comps(3);
end
