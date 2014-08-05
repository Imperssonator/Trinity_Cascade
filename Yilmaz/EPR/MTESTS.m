function out = MTESTS(mix)
% 100, 50 , 30, 10, 5wt% P3HT to start with 
% 1  , 2  , 3 , 4 , 5    are the indices of the M matrix with these comps
% you want M(mix,3) for P3HT VF, M(mix,4) for PS VF to start

steps = 10;
VFOs = zeros(3,steps);

%M = Mincheol();

% startP3HT = M(mix,3);
% startPS = M(mix,4);
startP3HT = 0.005; % for Prausnitz
startPS = 0.005; % for Prausnitz
startsolv = 1-startP3HT-startPS;
endsolv = 0.01;
solvcomps = linspace(endsolv,startsolv,steps);

for i = 1:steps
    VFOs(3,i) = solvcomps(i);
    VFOs(1,i) = (1-solvcomps(i))*startP3HT/(startP3HT+startPS);
    VFOs(2,i) = 1-VFOs(1,i)-VFOs(3,i);
end

TP = [];
% DP = [119 1635 1]'; %Mincheol
DP = [1000 1000 1]'; %Prausnitz

figure
hold on

for i = 1:steps
    disp(i)
    disp(VFOs(:,i))
    [VFeq,PFeq,Eeq,stab] = EPR(VFOs(:,i),DP,0);
    disp(VFeq)
    disp(PFeq)
    x = VFOs(3,i);
    if stab == 0
        plot(x,VFeq(1,1),'or',x,VFeq(1,2),'xr',x,VFeq(2,1),'ok',x,VFeq(2,2),'xk','MarkerSize',12)
        TP = [TP,VFeq];
    elseif stab == 1
        plot(x,VFeq(1,1),'ob',x,VFeq(1,2),'xb',x,VFeq(2,1),'oc',x,VFeq(2,2),'xc','MarkerSize',12)
        TP = [TP,VFeq];
    else
        plot(x,VFeq(1,1),'og',x,VFeq(2,1),'oy','MarkerSize',12)
    end
end

hold off
out = TP;
figure
ternplot(TP(1,:),TP(3,:),TP(2,:),'ob')

end