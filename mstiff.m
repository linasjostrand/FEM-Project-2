function D = mstiff(ep,eff)
%MSTIFF calculate material stiffness

E = ep(1); v = ep(2);
F = [eff(1:2)' 0; eff(3:4)' 0; 0 0 1];
K = E/(3*(1-2*v)); G = E/(2*(1+v));
J = det(F);
C = F'*F;
invC = inv(C);

a1 = K*J^2 + 2*G/9*J^(-2/3)*trace(C);
a2 = 2*G/3*J^(-2/3);
a3 = G/3*J^(-2/3)*trace(C)-K/2*(J^2-1);

D = zeros(3,3);
invC = [invC(1,1) invC(2,2) invC(1,2)];
% ij = 11 22 12
for ij = 1:3
    for kl = 1:3
        D(ij,kl) = a1*invC(ij)*invC(kl) - a2*((ij~=3)*invC(kl) + (kl~=3)*invC(ij));
    end
end

D(1,1) = D(1,1) + 2*a3*invC(1)^2;
D(1,2) = D(1,2) + 2*a3*invC(3)^2;
D(1,3) = D(1,3) + 2*a3*invC(1)*invC(3);
D(2,1) = D(2,1) + 2*a3*invC(3)^2;
D(2,2) = D(2,2) + 2*a3*invC(2)^2;
D(2,3) = D(2,3) + 2*a3*invC(3)*invC(2);
D(3,1) = D(3,1) + 2*a3*invC(1)*invC(3);
D(3,2) = D(3,2) + 2*a3*invC(3)*invC(2);
D(3,3) = D(3,3) + a3*(invC(1)*invC(2) + invC(3)*invC(3));


end