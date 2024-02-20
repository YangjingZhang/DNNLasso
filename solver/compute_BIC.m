function [llk,df] = compute_BIC(Q,P,Gamma,Omega,G,Theta)
p1 = size(Q,1);
p2 = size(P,2);
tmp1 = ones(p1,p1);
tmp1 = triu(tmp1,1);
ind1 = (tmp1(:) > 0);
tmp2 = ones(p2,p2);
tmp2 = triu(tmp2,1);
ind2 = (tmp2(:) > 0);
Omega = (Omega + Omega')/2;
Gamma = (Gamma + Gamma')/2;
Omega_d = eig(full(Omega));
Gamma_d = eig(full(Gamma));
if min(Omega_d) + min(Gamma_d) < 0
    fprintf('[not psd]');
end
llk = -sum(sum(log(Omega_d+Gamma_d'))) + sum(sum(Omega.*P)) + sum(sum(Gamma.*Q));
df = nnz(G(ind1)) + nnz(Theta(ind2));
end
