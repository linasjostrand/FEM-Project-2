function es = stresscal(ep, eff)
%STRESSCAL 
% ep = [E v]
% F deformation gradient

F = [eff(1:2)' 0; eff(3:4)' 0; 0 0 1];
E = ep(1); v = ep(2);
K = E/(3*(1-2*v)); G = E/(2*(1+v));
J = det(F);
C = F'*F;

S = K/2*(J^2-1)*eye(3)/C + G*J^(-2/3)*(eye(3) - 1/3*trace(C)*eye(3)/C);
es = [S(1,1); S(2,2); S(1,2)];

end