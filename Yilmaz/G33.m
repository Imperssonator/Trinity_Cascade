function out = G33(VF,DP,G)

u1 = VF(1)/(VF(1)+VF(2));
u2 = VF(2)/(VF(1)+VF(2));
dg12_du2 = 0;
dg13_dVF3 = 0;
dg23_dVF3 = 0;
d2g12_d2u2=0;
d2g13_d2VF3=0;
d2g23_d2VF3=0;

out = 1/VF(1)+DP(1)/(DP(3)*VF(3))-2*G(1,3)-2*u2^2*(1-u1)*dg12_du2+u1*u2^3*d2g12_d2u2...
    + 2*(VF(1)-VF(3))*dg13_dVF3+VF(1)*VF(3)*d2g13_d2VF3+2*DP(1)/DP(2)*VF(2)*dg23_dVF3...
    + DP(1)/DP(2)*VF(2)*VF(3)*d2g23_d2VF3;

end