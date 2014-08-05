function out = DTESTS()

steps = 5;
VFOs = zeros(3,steps);
j = 1;

for i = [5 10 15 20 30 40]
    [VFO,DP] = Dalsu(i,20);
    VFOs(:,j) = VFO;
    j = j+1;
end

TP = [];

figure
hold on

for i = 1:6
    disp(i)
    disp(VFOs(:,i))
    [VFeq,PFeq,Eeq,stab] = EPR(VFOs(:,i),DP,0);
    disp(VFeq)
    disp(PFeq)
    x = i;
    if stab == 0
        plot(x,VFeq(3,1),'or',x,VFeq(3,2),'xr','MarkerSize',12)
        TP = [TP,VFeq];
    elseif stab == 1
        plot(x,VFeq(3,1),'ob',x,VFeq(3,2),'xb','MarkerSize',12)
        TP = [TP,VFeq];
    else
        plot(x,VFeq(3,1),'og','MarkerSize',12)
    end
end

hold off
out = TP;
figure
ternplot(TP(1,:),TP(3,:),TP(2,:),'ob')

end