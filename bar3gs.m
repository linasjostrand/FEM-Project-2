function [es,ee] = bar3gs(ec,ep,ed)
%BAR3GS Compute the strain and normal force in a 3D bar element
%   [es,ee] = bar3gs(ec,ep,ed)

EA = prod(ep);
x0 = ec(:,2)-ec(:,1);
l0 = norm(x0);
u = ed(4:end)' - ed(1:3)';
x = x0 + u;
l = norm(x);
ee = (l^2-l0^2)/(2*l0^2);
es = EA*ee;




end