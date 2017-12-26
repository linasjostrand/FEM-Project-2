function N = norfb(ec,ee,k,r)
%NORFB Calculates normal force with special constitutive model

x0 = ec(:,2)-ec(:,1);
l0 = norm(x0);
lambda = sqrt(2*ee + 1);
lambda_c = r/l0;

N = min(k*(1-lambda_c/lambda), 0);


end