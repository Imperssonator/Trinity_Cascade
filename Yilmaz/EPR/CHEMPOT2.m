function out = CHEMPOT2(VF,DP,G)

u1 = VF(1)/(VF(1)+VF(2));
u2 = VF(2)/(VF(1)+VF(2));
g12_u2 = 0;
g13_VF3 = 0;
g23_VF3 = 0;

out = log(VF(2))+1-VF(2)-DP(2)/DP(1)*VF(1)-DP(2)/DP(3)*VF(3)+(G(1,2)*DP(1)/DP(2)*VF(1)...
    +G(2,3)*VF(3))*(VF(1)+VF(3))-G(1,3)*DP(2)/DP(1)*VF(1)*VF(3)+u1*u2*DP(2)/DP(1)*VF(1)*g12_u2...
    -DP(2)/DP(1)*VF(1)*VF(3)^2*g13_VF3-VF(2)*VF(3)^2*g23_VF3;

end