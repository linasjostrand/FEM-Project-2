function ef = bar3gf(ec,ed,es)
%BAR3GF Compute the internal element force vector for a 3D bar element
%   ef = bar3gf(ec,ed,es)

x0 = ec(:,2)-ec(:,1);
l0 = norm(x0);
u = ed(4:end)' - ed(1:3)';
x = x0 + u;
fA = -es/l0*x;
fB = es/l0*x;
ef = [fA; fB];

end

