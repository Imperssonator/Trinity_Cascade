VFOs = zeros(3,10);

startNS = 0.8;
endNS = 0.1;
steps = 10;
NScomps = linspace(endNS,startNS,steps);

for i = 1:steps
    VFOs(1,i) = NScomps(i);
    VFOs(2,i) = (1-VFOs(1,i))/2;
    VFOs(3,i) = 1-VFOs(1,i)-VFOs(2,i);
end

TP = zeros(1,3);
PFi = zeros(3,2)+0.5;
DP = [1 4 500]';

for i = 1:steps
    disp(i)
    [VFeq,PFeq,Eeq,stab] = EPR(VFOs(:,i),DP,0);
    if stab == 0 || stab == 1
        TP = [TP;VFeq'];
    end
end

ternplot(TP(:,1),TP(:,3),TP(:,2),'or')