function Y = prox_L1_off_diag_NN(X,lambda)
tmp = proj_inf_off(X,lambda);
Y = X - tmp;
Y = Y + diag(max(0,diag(X)));
end
function Y = proj_inf_off(X,lambda)
if lambda < 0 
   error('Lambda needs to be nonegative');
elseif lambda == 0
   Y = 0*X; 
else
   Y = max(-lambda,min(X,lambda));
end
Y = Y + diag(diag(X) - diag(Y));
end
