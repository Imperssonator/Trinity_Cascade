function Energy = Gibbs(N,VF,g,T)
%% GIBBS
% Finally, the real deal, a general formula for the Gibbs energy of a phase
% of arbitrary number of components with specified numbers of moles [N],
% volume fractions [VF], binary F-H interaction parameters in a 3x3 [g],
% and temperature T

Entropic = sum(N.*(log(VF)));
Enthalpic = 0;
components = length(N);

for i = 1:(components-1)
    for j = (i+1):components
        Enthalpic = Enthalpic + g(i,j)*N(i)*VF(j);
    end
end

Energy = 8.314*T*(Entropic+Enthalpic);

end
