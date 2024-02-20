function s = findzero(z,beta,mu,s0)
f = @(x) (x-z)/beta - sum(1./(x+mu));
epsilon = max(1/beta,1)*1e-12;
lft = -min(mu);
if nargin > 3
    s = s0;
else
    s = lft + 1;
end
maxiter = 10;
for i = 1:maxiter
    h = f(s);
    if abs(h) < epsilon
        break;
    end
    hdiff = 1/beta + sum(1./(s+mu).^2);
    s = s - h/hdiff; 
    if isnan(s) 
        break;
    end
end
if isnan(s) || (s <= -min(mu)) || (i == maxiter)
    if nargin > 3
        s = findzero_bisection(z,beta,mu,s0);
    else
        s = findzero_bisection(z,beta,mu);
    end
end
end
