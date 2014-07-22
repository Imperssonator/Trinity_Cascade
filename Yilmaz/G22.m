function out = G22(VF,DP,G)

u1 = VF(1)/(VF(1)+VF(2));
u2 = VF(2)/(VF(1)+VF(2));
dg12_du2 = 0;
d2g12_d2u2=0;

out = 1/VF(1)+DP(1)/(DP(2)*VF(2))-2*G(1,2)-2*(u1-u2)*dg12_du2+u1*u2*d2g12_d2u2;

end