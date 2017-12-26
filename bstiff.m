function D = bstiff(ec,ee,k,r)
%BSTIFF calculate material stiffness in bars
x0 = ec(:,2)-ec(:,1);
l0 = norm(x0);
lambda = sqrt(2*ee + 1); % stretch
lambda_c = r/l0;

D = (lambda < lambda_c)*k*lambda_c/(2*ee+1)^(3/2);

end