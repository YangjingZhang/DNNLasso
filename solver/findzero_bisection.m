function root = findzero_bisection(z,beta,mu,s0)
f = @(x) (x-z)/beta - sum(1./(x+mu));
epsilon = max(1/beta,1)*1e-12;

mu1 =  min(mu);
if nargin < 4
    lft = -mu1;
else
    lft = max(s0,-mu1);
end
dn = 0.1;
rht = lft + dn;
if f(rht) < 0
    while f(rht) < 0
        dn = dn*2;
        rht = rht + dn;
    end
    lft = rht - dn;
end

if abs(f(rht))< epsilon
    root = rht;
    return;
elseif abs(f(lft)-0.0)< epsilon
    root = lft;
    return;
end

root = subfun(f,lft,rht,epsilon);
end


function y = subfun(f,x0,x1,epsilon)
x2 = (x0 + x1) / 2;
err = abs(x0 - x1);
y = x2;
count = 0;
while err > epsilon 
    if f(x2) > 0
        x1 = x2;
    else
        x0 = x2;
    end
    x2 = (x0 + x1) / 2;
    err = min(abs(x2-x1),abs(x2-x0));%abs(x2-x1);
    y = x2;
    count = count + 1;
    if count > 10000; keyboard; end
end
end

