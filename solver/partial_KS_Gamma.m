function AtH = partial_KS_Gamma(Gamma_V,Gamma_d,Omega_V,Omega_d)
tmp = (sum(1./(Gamma_d + Omega_d'),2)).^0.5;
half_AtH = Gamma_V.*tmp';
AtH = half_AtH*half_AtH';
end


% function AtH = partial_KS_Gamma(Gamma_V,Gamma_d,Omega_V,Omega_d)
% p1 = length(Gamma_d);
% AtH = zeros(p1,p1);
% for i = 1:p1
%     AtH = AtH + sum(sum((Omega_V.*Omega_V).*(1./(Gamma_d(i)+Omega_d'))))*(Gamma_V(:,i)*Gamma_V(:,i)');
% end
% end