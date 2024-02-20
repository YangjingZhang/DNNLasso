%%*************************************************************************
%% DNNLasso:
%% Alternating direction method of multipliers for solving 
%% diagonally non-negative graphical lasso with Kronecker-sum precison matrix:
%% min -log|Omega + Gamma| + <Omega,P> + <Gamma,Q> + lamt*||Gamma||_1 + lams*||Omega||_1
%% where log| | is the log-determinant function
%%       || ||_1 sums the off-diagonal entries of a matrix
%%       lamt, lams are positive penalty parameters
%%       sample covariance Q \in R^{p1*p1}, P \in R^{p2*p2}
%%
%%*************************************************************************
%% DNNLasso: 
%% Copyright (c) 2024 by
%% Meixia Lin and Yangjing Zhang
%%*************************************************************************
%% Input - sample covariance: Q \in R^{p1*p1}, P \in R^{p2*p2}
%%       - penalty parameters: lamt>0, lams>0
%%       - OPTIONS.tol: tolerance, by default = 1e-6
%%       - OPTIOMS.maxiter: maximum iterations, by default = 10000
%%       - initilization: Gamma0,Omega0,G0,Theta0,X0,Y0, by default = (I_p1,I_p2,I_p1,I_p2,0,0)
%% Output - obj: optimal objective function value
%%        - Gamma, Omega: estimated precision matrices
%%        - G,Theta,X,Y: dual variables
%%        - info: additional information
%%        - runhist: running history
%%
function [obj,Gamma,Omega,G,Theta,X,Y,info,runhist] = DNNLasso(Q,P,lamt,lams,OPTIONS,Gamma0,Omega0,G0,Theta0,X0,Y0)
%%
if nargin < 5
    OPTIONS = [];
end
maxiter = 10000;
tol = 1e-6;
sig0 = 1;
sig_fix = 0;
printlevel = 1;
tau = 1.618;
info = [];
prim_flag = 1; 
dual_flag = 1;
if isfield(OPTIONS,'tol');          tol           = OPTIONS.tol;        end
if isfield(OPTIONS,'printlevel');   printlevel    = OPTIONS.printlevel; end
if isfield(OPTIONS,'maxiter');      maxiter      = OPTIONS.maxiter;     end
p1 = size(Q,1);
p2 = size(P,1);
%print dimensions
if printlevel
    fprintf('\n dim. of temporal sdp var = %2.0d',p1);
    fprintf('\n dim. of spatial sdp var = %2.0d',p2);
    fprintf('\n lamt = %2.1e, lams = %2.1e\n',lamt,lams);
end
tstart = clock;
normQ = norm(Q,'fro');
normP = norm(P,'fro');
%% initialization
if ~exist('Gamma0','var') || ~exist('Omega0','var') || ~exist('G0','var') || ~exist('Theta0','var') || ~exist('X0','var') || ~exist('Y0','var')
    Gamma = eye(p1);
    Omega = eye(p2);
    G = Gamma;
    Theta = Omega;
    Omega_d = ones(p2,1);
    X = zeros(p1,p1);
    Y = zeros(p2,p2);
else
    Gamma = Gamma0;
    Omega = Omega0;
    G = G0;
    Theta = Theta0;
    [~,Omega_D] = mexeig(full(Omega));
    Omega_d = diag(Omega_D);
    X = X0;
    Y = Y0;
end
msg = '';
breakyes = 0;
prim_win = 0;
dual_win = 0;
new_flag = 0;
fprintf('\n*****************************************************')
fprintf('*************************************')
fprintf('\n  DNNLasso')
fprintf('\n*****************************************************')
fprintf('*************************************')
fprintf('\n  iter   |  pinfeas    dinfeas   compliment  | relgaporg |      pobj            dobj      |')
fprintf('  time  |  sigma   |')
%% main loop
for iter = 1:maxiter
    %update Gamma
    Gammatmp = G + (1/sig0)*(X - Q);
    [Gamma_V,Gammatmp_D] = eig(full(Gammatmp)); 
    Gammatmp_d = diag(Gammatmp_D);
    Gamma_d = zeros(p1,1);
    Gamma_d(1) = findzero(Gammatmp_d(1),1/sig0,Omega_d);
    for i = 2:p1
        Gamma_d(i) = findzero(Gammatmp_d(i),1/sig0,Omega_d,Gamma_d(i-1));
    end
    tmp = Gamma_d.*Gamma_V';
    Gamma = tmp'*Gamma_V';
    Gamma = (Gamma + Gamma')/2;
    %update Omega
    Omegatmp = Theta + (1/sig0)*(Y - P);
    [Omega_V,Omegatmp_D] = eig(full(Omegatmp)); 
    Omegatmp_d = diag(Omegatmp_D);
    Omega_d = zeros(p2,1);
    Omega_d(1) = findzero(Omegatmp_d(1),1/sig0,Gamma_d);
    for i = 2:p2
        Omega_d(i) = findzero(Omegatmp_d(i),1/sig0,Gamma_d,Omega_d(i-1));
    end
    tmp = Omega_d.*Omega_V';
    Omega = tmp'*Omega_V';
    Omega = (Omega + Omega')/2;
    %update G, Theta
    G = prox_L1_off_diag_NN(Gamma-X/sig0,lamt/sig0);
    Theta = prox_L1_off_diag_NN(Omega-Y/sig0,lams/sig0);
    %update X, Y
    X = X - (tau*sig0)*(Gamma - G);
    Y = Y - (tau*sig0)*(Omega - Theta);
    %Compute kkt residual
    normGamma = norm(Gamma,'fro');
    normOmega = norm(Omega,'fro');
    normG = norm(G,'fro');
    normTheta = norm(Theta,'fro');
    Rp1 = norm(Gamma-G,'fro')/(1+normGamma+normG);
    Rp2 = norm(Omega-Theta,'fro')/(1+normOmega+normTheta);
    prim_infeas = max(Rp1,Rp2);
    partial_Gamma = partial_KS_Gamma(Gamma_V,Gamma_d,Omega_V,Omega_d);
    partial_Omega = partial_KS_Omega(Gamma_V,Gamma_d,Omega_V,Omega_d);
    Rd1 = norm(-partial_Gamma+Q-X,'fro')/(1+norm(partial_Gamma,'fro')+normQ+norm(X,'fro'));
    Rd2 = norm(-partial_Omega+P-Y,'fro')/(1+norm(partial_Omega,'fro')+normP+norm(Y,'fro'));
    dual_infeas = max(Rd1,Rd2);
    proxGX = prox_L1_off_diag_NN(G-X,lamt);
    com1 = norm(proxGX-G,'fro')/(1+norm(proxGX,'fro')+normG);
    proxThetaY = prox_L1_off_diag_NN(Theta-Y,lams);
    com2 = norm(proxThetaY-Theta,'fro')/(1+norm(proxThetaY,'fro')+normTheta);
    complimentarity = max(com1,com2);
    maxinfeas = max([prim_infeas,dual_infeas,complimentarity]);
    %compute the objective function
    tmp = sum(sum(log(Omega_d+Gamma_d')));
    if iter > 1
        oldobj = obj(1);
    else
        oldobj = 0;
    end
    obj(1) = - tmp + sum(sum(Omega.*P)) + sum(sum(Gamma.*Q)) + lamt*(sum(sum(abs(Gamma)))-sum(abs(diag(Gamma)))) + lams*(sum(sum(abs(Omega)))-sum(abs(diag(Omega))));
    obj(2) = - tmp + sum(sum(Omega.*(P-Y))) + sum(sum(Gamma.*(Q-X)));
    %record running history
    gap = obj(1) - obj(2);
    rel_gap = abs(gap)/(1+sum(abs(obj)));
    ttime = etime(clock,tstart);
    runhist.prim_obj(iter)    = obj(1);
    runhist.dual_obj(iter)    = obj(2);
    runhist.obj_ratio(iter) = abs((obj(1)-oldobj)/obj(1));
    runhist.gap(iter)         = gap;
    runhist.relgap(iter)      = rel_gap;
    runhist.prim_infeas(iter) = prim_infeas;
    runhist.dual_infeas(iter) = dual_infeas;
    runhist.feasratio(iter) = prim_infeas/dual_infeas;
    runhist.complimentarity(iter) = complimentarity;
    runhist.maxinfeas(iter)   = maxinfeas;
    runhist.cputime(iter) = ttime;
    runhist.sigma(iter)   = sig0;
    %check for termination
    if maxinfeas < tol
        msg = sprintf('max(prim_infeas,dual_infeas,complimentarity) < %3.2e',tol);
        breakyes = 1;
    end
    %print results
    if printlevel  && (iter <= 20 || rem(iter,100) == 1 || breakyes || iter == maxiter)
        fprintf('\n   %4.0d  |  %3.2e   %3.2e   %3.2e   |  %3.2e | %- 8.7e %- 8.7e  |  %3.2f  | %3.2e |',...
            iter,prim_infeas,dual_infeas,complimentarity,rel_gap,obj(1),obj(2),ttime,sig0);
        if iter > 1
            fprintf('relobj = %3.2e', (obj(1) - oldobj)/obj(1));
        end
    end
    if breakyes; break; end
    %update sigma
    if runhist.feasratio(iter) < 1
        prim_win = prim_win+1;
    elseif runhist.feasratio(iter) >= 1
        dual_win = dual_win+1;
    end
    sigma_update_iter = sigma_fun(iter);
    sigmascale = 1.5;%3;%1.5;%1.25;
    if ~new_flag
        if (~sig_fix) && (rem(iter,sigma_update_iter)==0)
            sigmamax = 1e4; sigmamin = 1e-4;
            if (iter <= 4*2500)
                if (prim_win > max(1,1.2*dual_win))
                    prim_win = 0;
                    sig0 = max(sigmamin,sig0/3);
                elseif (dual_win > max(1,1.2*prim_win))
                    dual_win = 0;
                    sig0 = min(sigmamax,sig0*1.25);
                end
            else
                feasratiosub = runhist.feasratio(max(1,iter-19):iter);
                meanfeasratiosub = mean(feasratiosub);
                if meanfeasratiosub < 0.1 || meanfeasratiosub > 1/0.1
                    sigmascale = 1.4;
                elseif meanfeasratiosub < 0.2 || meanfeasratiosub > 1/0.2
                    sigmascale = 1.35;
                elseif meanfeasratiosub < 0.3 || meanfeasratiosub > 1/0.3
                    sigmascale = 1.32;
                elseif meanfeasratiosub < 0.4 || meanfeasratiosub > 1/0.4
                    sigmascale = 1.28;
                elseif meanfeasratiosub < 0.5 || meanfeasratiosub > 1/0.5
                    sigmascale = 1.26;
                end
                primidx = find(feasratiosub <= 1);
                dualidx = find(feasratiosub >  1);
                if (length(primidx) >= 12)
                    if prim_flag
                        sig0 = max(sigmamin,sig0/sigmascale);
                    end
                end
                if (length(dualidx) >= 12)
                    if dual_flag
                        sig0 = min(sigmamax,sig0*sigmascale);
                    end
                end
            end
        end
    else
        if (~sig_fix) && (rem(iter,sigma_update_iter)==0)
            sigmamax = 1e4; sigmamin = 1e-4;
            if (iter <= 4*2500)
                if (prim_win > max(1,1.2*dual_win))
                    if prim_flag
                        prim_win = 0;
                        sig0 = max(sigmamin,sig0/3);
                        dual_flag = 0;
                    end
                elseif (dual_win > max(1,1.2*prim_win))
                    if dual_flag
                        dual_win = 0;
                        sig0 = min(sigmamax,sig0*1.25);
                        prim_flag = 0;
                    end
                end
            else
                feasratiosub = runhist.feasratio(max(1,iter-19):iter);
                meanfeasratiosub = mean(feasratiosub);
                if meanfeasratiosub < 0.1 || meanfeasratiosub > 1/0.1
                    sigmascale = 1.4;
                elseif meanfeasratiosub < 0.2 || meanfeasratiosub > 1/0.2
                    sigmascale = 1.35;
                elseif meanfeasratiosub < 0.3 || meanfeasratiosub > 1/0.3
                    sigmascale = 1.32;
                elseif meanfeasratiosub < 0.4 || meanfeasratiosub > 1/0.4
                    sigmascale = 1.28;
                elseif meanfeasratiosub < 0.5 || meanfeasratiosub > 1/0.5
                    sigmascale = 1.26;
                end
                primidx = find(feasratiosub <= 1);
                dualidx = find(feasratiosub >  1);
                if (length(primidx) >= 12)
                    if prim_flag
                        sig0 = max(sigmamin,sig0/sigmascale);
                        dual_flag = 0;
                    end
                end
                if (length(dualidx) >= 12)
                    if dual_flag
                        sig0 = min(sigmamax,sig0*sigmascale);
                        prim_flag = 0;
                    end
                end
            end
        end
    end

    if ~new_flag
        objratiosub = runhist.obj_ratio(max(1,iter-19):iter);
        meanobjratiosub = mean(objratiosub);
        if meanobjratiosub < 1e-9
            new_flag = 1;
            if sig0 < 1
                prim_flag = 1; dual_flag = 0; %dual_flag = 1; prim_flag = 0;
            else
                dual_flag = 1; prim_flag = 0; %prim_flag = 1; dual_flag = 0;
            end
            sig0 = max(1,sig0);
            prim_win = 0;
            dual_win = 0;
            sig_fix = 1;
        end
    end
end %end of main loop
%%
info.iter    = iter;
info.obj     = [obj(1),obj(2)];
info.gap     = abs(obj(1)-obj(2));
info.relgap  = rel_gap;
info.pinfeas = prim_infeas;
info.dinfeas = dual_infeas;
info.complimentarity = complimentarity;
info.time    = ttime;
info.msg     = msg;
if (maxinfeas < tol)
    info.termcode = 0;
elseif (prim_infeas < sqrt(tol) && dual_infeas < tol)
    msg = sprintf('problem is partially solved to the required tolerance');
    info.termcode = -1;
else
    if (iter == maxiter)
        msg = sprintf('maximum number of iterations reached');
        info.termcode = 2;
    end
end
fprintf('\n--------------------------------------------------------')
fprintf('\n %s',msg);
fprintf('\n--------------------------------------------------------')
fprintf('\n primal objval = %9.8e',obj(1));
fprintf('\n dual   objval = %9.8e',obj(2));
fprintf('\n gap = %3.2e',gap);
fprintf('\n relative gap = %3.2e',rel_gap);
fprintf('\n prim_infeas = %3.2e',prim_infeas);
fprintf('\n dual_infeas = %3.2e',dual_infeas);
fprintf('\n complimentarity = %3.2e',complimentarity);
ttime = etime(clock,tstart);
fprintf('\n time = %3.2f',ttime);
fprintf('\n--------------------------------------------------------')
fprintf('------------------------')
fprintf('\n')
end
%%
function sigma_update_iter = sigma_fun(iter)
if (iter < 30)
    sigma_update_iter = 10;%6;%3;
elseif (iter < 60)
    sigma_update_iter = 12;%6;
elseif (iter < 120)
    sigma_update_iter = 12;
elseif (iter < 250)
    sigma_update_iter = 25;
elseif (iter < 500)
    sigma_update_iter = 50;
elseif (iter < 2000)
    sigma_update_iter = 100;%200;
elseif (iter < inf)  
    sigma_update_iter = 500;
end
end
