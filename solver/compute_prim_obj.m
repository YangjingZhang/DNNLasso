function [fx,flag] = compute_prim_obj(Gamma,Omega,P,Q,lamt,lams)
Omega = (Omega + Omega')/2;
Gamma = (Gamma + Gamma')/2;
Omega_d = eig(full(Omega));
Gamma_d = eig(full(Gamma));
flag = 1;
if min(Omega_d)+min(Gamma_d) < 0
    fprintf('[not psd]');
    flag = 0;
%     keyboard;
end
tmp = sum(sum(log(Omega_d+Gamma_d'))) ;
fx = - tmp + sum(sum(Omega.*P)) + sum(sum(Gamma.*Q)) + lamt*(sum(sum(abs(Gamma)))-sum(abs(diag(Gamma)))) + lams*(sum(sum(abs(Omega)))-sum(abs(diag(Omega))));
end