function ef = plan3gf(ec,t,ed,es)
%PLAN3GF computes the internal element forces vector

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

ef = B0'*es*A*t;




end