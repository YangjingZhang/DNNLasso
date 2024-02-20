function BtH = partial_KS_Omega(Gamma_V,Gamma_d,Omega_V,Omega_d)
tmp = (sum(1./(Omega_d + Gamma_d'),2)).^0.5;
half_BtH = Omega_V.*tmp';
BtH = half_BtH*half_BtH';
end


% function BtH = partial_KS_Omega(Gamma_V,Gamma_d,Omega_V,Omega_d)
% p1 = length(Gamma_d);
% p2 = length(Omega_d);
% BtH = zeros(p2,p2);
% for i = 1:p1
%     BtH = BtH +  (Omega_V.*(1./(Gamma_d(i)+Omega_d')))*Omega_V';
% end
% end