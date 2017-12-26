function [ee,eff] = plan3gs(ec,ed)
%PLAN3GS computes the Green-Lagrange strains and the deformation gradient
% ee strains
% eff deformation gradient
% ec element nodal coordinates in undeformed config
% current element displacement vector

ex = ec(1,:); ey = ec(2,:);
C=[ 1  ex(1) ey(1)   0     0       0
    0    0     0     1   ex(1)   ey(1)
    1  ex(2) ey(2)   0     0       0
    0    0     0     1   ex(2)   ey(2)
    1  ex(3) ey(3)   0     0       0
    0    0     0     1   ex(3)   ey(3)];
% Nonlinear part of B0 matrix
alp = C\ed;
dudx = [0 1 0 0 0 0
        0 0 0 0 1 0]*alp;
dudy = [0 0 1 0 0 0
        0 0 0 0 0 1]*alp;

% Displacement gradient tensor
D = [dudx,dudy];
% Deformation gradient
F = D + eye(2);
eff = [F(1,:) F(2,:)]'; %vector format

% Green Lagrange strains
E = 0.5*(F'*F-eye(2));
ee = [E(1,1) E(2,2) 2*E(1,2)]';


end