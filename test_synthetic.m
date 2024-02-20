%% On synthetic type1 or type2 graphs
%% plot obj vs time
HOME = pwd;
addpath(genpath(HOME));
close all; clear; rng('default');
p1 = 500; p2 = 100; %problem dimensions
type = [2,2]; %graph type
n = 1; %sample size
lams = p1*0.01;%penalty parameters
lamt = p2*0.01;
Omega0 = generate_precision(p2,type(1));
Gamma0 = generate_precision(p1,type(2));
Omega0_off = Omega0 - diag(diag(Omega0));
Gamma0_off = Gamma0 - diag(diag(Gamma0));
[U_Omega,D_Omega] = eig(Omega0);
[U_Gamma,D_Gamma] = eig(Gamma0);
D = 1./(diag(D_Gamma)+diag(D_Omega)');
Dhalf = D.^0.5;
%generate sample covariance matrices Q and P
Q = zeros(p1,p1);
P = zeros(p2,p2);
for k = 1:n
    Rtmp = Dhalf.*randn(p1,p2);
    Z = U_Gamma*Rtmp*U_Omega';
    Q = Q + Z*Z';
    P = P + Z'*Z;
end
Q = Q/n;
P = P/n;
%% benchmark: run DNNLasso with tol=1e-8 to obtain the 'true' obj value
OPTIONS.tol = 1e-8;
OPTIONS.maxiter = 10000;
[~,Gamma,Omega,~,~,~,~,~,runhist_admm] = DNNLasso(Q,P,lamt,lams,OPTIONS);
optimal_obj = min(runhist_admm.prim_obj);
record0 = zeros(8,10000);
record0(8,5) = optimal_obj;
%% run DNNLasso
OPTIONS.tol = 1e-7; %tolerance
OPTIONS.maxiter = 10000; %max iterations
[obj,~,~,Gamma_admm,Omega_admm,~,~,info_admm,runhist_admm] = DNNLasso(Q,P,lamt,lams,OPTIONS);
iter_admm = info_admm.iter;
record0(8,3) = iter_admm;
record0(1,1:iter_admm) = runhist_admm.prim_obj(1:iter_admm);
record0(2,1:iter_admm) = runhist_admm.cputime(1:iter_admm);
record0(3,1:iter_admm) = runhist_admm.maxinfeas(1:iter_admm);
%% plot obj vs time
figure;
iter_admm = record0(8,3);
optimal_obj = record0(8,5);
obj_admm = record0(1,1:iter_admm);
time_admm = record0(2,1:iter_admm);
kkt_admm = record0(3,1:iter_admm);
semilogy(time_admm(1:iter_admm),(obj_admm(1:iter_admm) - optimal_obj)/abs(optimal_obj),'LineWidth',1,'color','r');
ftsize = 15;
xlabel('time (sec)','fontsize',ftsize);
ylabel('relative objective','fontsize',ftsize);
axis square;
legend('DNNLasso','fontsize',ftsize,'location','northeast');
set(gca,'FontSize',ftsize);
set(gcf,'Position',[50 50 550 550]);
