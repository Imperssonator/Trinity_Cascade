function out = DGRAM()
% Make the binodal curve for CHCl3/hexane/P3HT by RANDOMLY TRYING POINTS

points = 50;
VFOs = zeros(3,points);

for i = 1:points
    r1 = rand; r2 = rand;
    VFOs(1,i) = r1;
    VFOs(2,i) = (1-r1)*r2;
    VFOs(3,i) = 1-VFOs(1,i)-VFOs(2,i);
end

TP = [];
DP = [1 1.6 50]'; % 39.1kD P3HT

for i = 1:points
    disp(i)
    disp(VFOs(:,i))
    [VFeq,PFeq,Eeq,stab] = EPR(VFOs(:,i),DP,0);
    disp(VFeq)
    disp(PFeq)
    if stab == 0
        %plot(x,VFeq(1,1),'or',x,VFeq(1,2),'xr',x,VFeq(2,1),'ok',x,VFeq(2,2),'xk','MarkerSize',12)
        TP = [TP,VFeq];
    elseif stab == 1
        %plot(x,VFeq(1,1),'ob',x,VFeq(1,2),'xb',x,VFeq(2,1),'oc',x,VFeq(2,2),'xc','MarkerSize',12)
        TP = [TP,VFeq];
    else
        %plot(x,VFeq(1,1),'og',x,VFeq(2,1),'oy','MarkerSize',12)
    end
end

out = TP;
figure
ternplot(TP(1,:),TP(3,:),TP(2,:),'ob')

end