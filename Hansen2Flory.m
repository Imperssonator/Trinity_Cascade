function out = Hansen2Flory(plot_temp)

%Hansen solubility parameters for P3HT, chloroform, Acetone, DCB, and THF
global ddP3HT dpP3HT dhbP3HT k R

ddP3HT = 19.45;
dpP3HT = 3.97; %MPa^1/2
dhbP3HT = 4.19;

R = 8.314;
k = 1.38E-23;

out = [XPChlo(298) XPAce(298)];
return

F_T = zeros(1000,2,100); %initialize free energy matrix, i = composition, k = Temp, j = 2 is ?F of mixing

for i = 1:1000
    F_T(i,1,:) = i/1000;
end

%polymer details
DOP = 240;

for i = 1:100 %this means T starts at 274 and goes to 374
    for j = 1:1000
        F_T(j,2,i) = DeltaF(F_T(j,1,i),DOP,i+273); %calculate ?F mix for a certain temperature/composition
    end
end

figure
plot(F_T(:,1,plot_temp-273),F_T(:,2,plot_temp-273))
out = F_T;

%disp(alpha*vChlo/R/T)

end

function out = DeltaF(c1,N,T)

global R
out = R*T*((c1/N)*log(c1)+(1-c1)*log(1-c1)+2*c1*(1-c1));

end

function out = XPChlo(T)

global ddP3HT dpP3HT dhbP3HT R

ddCHCl3 = 17.8;
dpCHCl3 = 3.1;
dhbCHCl3 = 5.7;
vChlo = 80.5; %cc/mol
alpha = 0.7;

out = alpha*vChlo/R/T*((ddP3HT-ddCHCl3)^2 + 0.25*(dpP3HT-dpCHCl3)^2 + 0.25*(dhbP3HT-dhbCHCl3)^2);

end

function out = XPAce(T)

global ddP3HT dpP3HT dhbP3HT R

ddAce = 15.50;
dpAce = 10.4;
dhbAce = 7;
vAce = 74;
alpha = 0.7;

out = alpha*vAce/R/T*((ddP3HT-ddAce)^2 + 0.25*(dpP3HT-dpAce)^2 + 0.25*(dhbP3HT-dhbAce)^2);

end

function out = XPDCB(T)

global ddP3HT dpP3HT dhbP3HT R

ddDCB = 19.2;
dpDCB = 6.3;
dhbDCB = 3.3;
vDCB = 112.8;
alpha = 1;

out = alpha*vDCB/R/T*((ddP3HT-ddDCB)^2 + 0.25*(dpP3HT-dpDCB)^2 + 0.25*(dhbP3HT-dhbDCB)^2);

if T==298
    %disp(out)
end

end

function out = XPTHF(T)

global ddP3HT dpP3HT dhbP3HT R

ddTHF = 16.8;
dpTHF = 5.7;
dhbTHF = 8.0;
vTHF = 81.7;
alpha = 1;

out = alpha*vTHF/R/T*((ddP3HT-ddTHF)^2 + 0.25*(dpP3HT-dpTHF)^2 + 0.25*(dhbP3HT-dhbTHF)^2);

end
