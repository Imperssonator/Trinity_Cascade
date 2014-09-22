function VFspin = BINSPIN(DP,VFO,VFeq,species,T)

npts = 1E5; %level of discretization of the volume fraction space
PTS = linspace(VFeq(1,1),VFeq(1,2),npts);
if PTS(2)-PTS(1) < 0
    PTS = fliplr(PTS);
end

GMAT = zeros(4,npts); %set up matrix to store polymer VFs, Gibbs energy values, second derivs
X = Chi(species,T);
GMAT(1,:) = PTS;
GMAT(2,:) = arrayfun(@(vf) MOLAR_GIBBS(vf,DP,X),GMAT(1,:));

h = 1/npts;
for ii = 2:npts-1
    GMAT(3,ii) = (GMAT(2,ii-1)+GMAT(2,ii+1)-2*GMAT(2,ii))/h^2;
end

VF1 = find(GMAT(3,:)<0,1,'first');
VF2 = find(GMAT(3,:)<0,1,'last');

VFspin = [PTS(VF1) PTS(VF2);...
          1-PTS(VF1) 1-PTS(VF2)];
end
