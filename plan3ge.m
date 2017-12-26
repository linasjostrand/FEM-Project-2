function Ke = plan3ge(ec,t,D,ed,es)
%PLAN3GE provides element stiffness matrix Ke for a triangular 3 node large
%deformation element in plane strain or stress.
%
% ec element nodal coordinates in undeformed config
% t element thickness
% ed current element displacement vector
% es 2d Piola-Kirchhoff stress tensor
% D constitutive matrix supplying material properties

ex = ec(1,:); ey = ec(2,:);
C=[ 1  ex(1) ey(1)   0     0       0  
    0    0     0     1   ex(1)   ey(1)
    1  ex(2) ey(2)   0     0       0  
    0    0     0     1   ex(2)   ey(2)
    1  ex(3) ey(3)   0     0       0  
    0    0     0     1   ex(3)   ey(3)];

A=1/2*det([ones(3,1) ex' ey']); % area
% Linear part of B0 matrix
B0l = [0 1 0 0 0 0
       0 0 0 0 0 1
       0 0 1 0 1 0]/C;

% Nonlinear part of B0 matrix
alp = C\ed;
Au = [alp(2) 0 alp(5) 0
      0 alp(3) 0 alp(6)
      alp(3) alp(2) alp(6) alp(5)];
H0 = [0 1 0 0 0 0
      0 0 1 0 0 0 
      0 0 0 0 1 0
      0 0 0 0 0 1]/C;

% Total B0 matrix
B0 = B0l+Au*H0;

% R matrix due to stresses
S = [es(1) es(3)
     es(3) es(2)];
R = blkdiag(S,S);


Ke=(B0'*D*B0 + H0'*R*H0)*A*t;




end