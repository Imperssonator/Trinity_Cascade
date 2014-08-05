function [X, HR] = Xij(solv,pol,T)

Hansen = struct('P3HT',[19.45 3.97 4.19 0]);
Hansen.CHCl3 = [17.8 3.1 5.7 80.5];
Hansen.Acetone = [15.5 10.4 7 74];
Hansen.DCB = [19.2 6.3 3.3 112.8];
Hansen.THF = [16.8 5.7 8.0 81.7];
Hansen.PS = [21.3 5.8 4.3 0];
Hansen.hexane = [14.9 0 0 131.6];

HV1 = Hansen.(solv);
HV2 = Hansen.(pol);
alpha = 1;
R = 8.314;

X = alpha*HV1(4)/R/T*((HV2(1)-HV1(1))^2 + 0.25*(HV1(2)-HV2(2))^2 + 0.25*(HV1(3)-HV2(3))^2);
HR = sqrt((HV2(1)-HV1(1))^2 + 0.25*(HV1(2)-HV2(2))^2 + 0.25*(HV1(3)-HV2(3))^2);

end