function Ke = bar3ge(ec,ep,ed,es)
%BAR3GE Compute stiffness matrix for 3D bar element
%   Compute stiffness matrix for 3D bar element. The element can be used 
%   for large deformations and rotations and is based on Green-Lagrange's
%   strain tensor.
%   Ke=bar3ge(ec,ep,ed,es)

EA = prod(ep);
x0 = ec(:,2)-ec(:,1);
l0 = norm(x0);
u = ed(4:end)' - ed(1:3)';
K0 = EA/l0^3*[x0*x0' -x0*x0'; -x0*x0' x0*x0'];
a = x0*u' + u*x0' + u*u';
Ku = EA/l0^3*[a -a; -a a];
Ks = es/l0*[eye(3) -eye(3); -eye(3) eye(3)];

Ke = K0 + Ku + Ks;

end