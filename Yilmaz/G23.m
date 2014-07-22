function out = G23(VF,DP,G)

u1 = VF(1)/(VF(1)+VF(2));
u2 = VF(2)/(VF(1)+VF(2));
dg12_du2 = 0;
dg13_dVF3 = 0;
dg23_dVF3 = 0;
d2g12_d2u2=0;

out = 1/VF(1)-(G(1,2)+G(1,3))+DP(1)/DP(2)*G(2,3)+u2*(u1-2*u2)*dg12_du2...
    + u1*u2^2*d2g12_d2u2-VF(3)*dg13_dVF3+DP(1)/DP(2)*VF(3)*dg23_dVF3;

end