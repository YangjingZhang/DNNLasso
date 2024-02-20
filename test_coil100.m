%% On coil100 data
%% object = 1: a rotating box of cold medicine
%% object = 100: a rotating box of a cargo
HOME = pwd;
addpath(genpath(HOME));
close all; clear;
object = 1;  %=1or100
scale = 1/16; %1/4, 1/16, reduce resolution
Z = zeros(72,(128*scale)^2);

    Q = zeros(72,72);
    P = zeros((128*scale)^2,(128*scale)^2);
    n = 3;
    for k = 1:n
        for pose = 0:71
            filename = ['coil100_',num2str(object),'_',num2str(pose*5),'.mat'];
            load(filename);
            coil100 = coil100(:,:,k);
            coil100 = imresize(coil100,scale);
            Z(pose+1,:) = coil100(:)';
        end
        Q = Q + Z*Z';
        P = P + Z'*Z;
    end
    Q = Q/n;
    P = P/n;

%% select penalty parameter by BIC
p1 = size(Q,1);
p2 = size(P,1);
rhobar_list = 10.^(-3.5:0.1:-2.0); %penalty parameters
record0 = zeros(2,length(rhobar_list));
select_again = 0;
if select_again
    for k = 1:length(rhobar_list)
        fprintf('\n %d-th of all %d', k, length(rhobar_list));
        lams = p1*rhobar_list(k);
        lamt = p2*rhobar_list(k);
        OPTIONS.tol = 5e-3;
        [obj,Gamma_admm,Omega_admm,G,Theta] = DNNLasso(Q,P,lamt,lams,OPTIONS);
        [record0(1,k),record0(2,k)] = compute_BIC(Q,P,Gamma_admm,Omega_admm,G,Theta);
    end
    figure(10); %plot Figure5(a)
    df_list1 = record0(2,:);
    nnzratio = df_list1/(p1*(p1-1)/2+p2*(p2-1)/2);
    tt = [1:length(rhobar_list)];
    llk_list1 = record0(1,:);
    gamma = 0.1;
    BIC_list1 = llk_list1 + (log(n)/n + 4*gamma*log(p1*p2))*df_list1;
    ftsize = 15;
    yyaxis left;
    semilogx(rhobar_list(tt),BIC_list1(tt),'-d','MarkerSize',5,'LineWidth',3);
    ylabel('BIC','fontsize',ftsize);
    yyaxis right;
    semilogx(rhobar_list(tt),nnzratio(tt)*100,':o','MarkerSize',5,'LineWidth',3);
    ytickformat('percentage');
    ylabel('sparsity level','fontsize',ftsize);
    xlabel('$\lambda_0$','fontsize',ftsize);
    xlim([rhobar_list(tt(1)),rhobar_list(tt(end))]);
    set(gca,'xtick',[]);
    xtick = rhobar_list(tt(1:2:end));
    xticklab = cellstr(num2str((-log10(xtick(:))), '10^-^{%1.1f}'));
    set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');
    legend('BIC','sparsity level');
    axis square;
    set(gca,'FontSize',ftsize);
    set(gcf,'Position',[0 0 500 500]);
    rl = rhobar_list(tt);
    bl = BIC_list1(tt);
    [~,dd] = min(bl);
    fprintf('\nBest ind is %d\n',dd);
    lam0_list = rl(dd);
else
    lam0_list = 10^(-2.7);
end
lam0 = lam0_list(1);
lamt = lam0*p2;
lams = lam0*p1;

    %%  DNNLasso
    OPTIONS.tol = 1e-6;
    tstart = clock;
    [obj,Gamma_admm,Omega_admm,G_admm,Theta_admm,X_admm,Y_admm,info_admm,runhist_admm] = DNNLasso(Q,P,lamt,lams,OPTIONS);
    iter_admm = info_admm.iter;
    Gamma = Gamma_admm - diag(diag(Gamma_admm));
    Omega = Omega_admm - diag(diag(Omega_admm));
    Theta_admm = Theta_admm - diag(diag(Theta_admm));


%% plot  obj vs time
iter_admm = info_admm.iter;
plot_seq_admm = [1:iter_admm];
figure;
plot(runhist_admm.cputime(plot_seq_admm),runhist_admm.prim_obj(plot_seq_admm),'LineWidth',3,'color','r');
ftsize = 15;
legend('DNNLasso','fontsize',ftsize,'location','northeast');
xlabel('time (sec)','fontsize',ftsize);
ylabel('objective','fontsize',ftsize);
axis square;
set(gca,'FontSize',ftsize);
set(gcf,'Position',[500 500 550 550]);
%% plot
figure(2);
G = G_admm;
G = abs(G);
G = G - diag(diag(G));
xx = 1:p1;
xx = xx - 0.5;
yy = 1:p1;
yy = yy - 0.5;
[X,Y] = meshgrid(xx,yy);
X = X(:);
Y = Y(:);
gg = G(:);
posind = gg > 0;
scatter(X(posind),Y(posind),gg(posind)*1.5,'ks','filled');
set(gca, 'YDir','reverse')
xlim([0 p1]);
ylim([0 p1]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
axis on;
box on;
axis square;
set(gcf,'Position',[500 0 550 550]);

%% plot relationship graph of frames from different angles
Tall = G_admm;
G0 = graph(Tall,'omitselfloops');
LWidths = abs(G0.Edges.Weight)/max(abs(G0.Edges.Weight));
LWidths = ones(size(LWidths));
name_pose = cell(72,1);
for pose = 0:71
    name_pose{pose+1,1} = sprintf('%d%c',pose*5, char(176));
end
G0.Nodes.Name = string(name_pose);
figure;
ppos = plot(G0,'LineWidth',LWidths);
layout(ppos,'circle');
ppos.MarkerSize = 5;
set(gca,'xticklabel','','yticklabel','');
axis square;
axis off;
axis tight;
ax=axis;
dxRange=(ax(2)-ax(1))/500;
dyRange=(ax(4)-ax(3))/500;
axis([ax(1)-dxRange,ax(2)+dxRange,ax(3)-dyRange,ax(4)+dyRange]);
rectangle('Position',[0.86 -0.05 0.3 0.1],'EdgeColor','r','LineWidth',3);
rectangle('Position',[-1.21 -0.05 0.3 0.1],'EdgeColor','r','LineWidth',3);
rectangle('Position',[0.48 0.81 0.12 0.3],'EdgeColor','r','LineWidth',3);
rectangle('Position',[-0.62 -1.15 0.12 0.3],'EdgeColor','r','LineWidth',3);
axes('pos',[.8 .45 .15 .15])
filename = ['coil100_',num2str(object),'_',num2str(0),'.mat'];
load(filename);
coil100 = imresize(coil100,scale);
imshow(coil100)
axes('pos',[.08 .45 .15 .15])
filename = ['coil100_',num2str(object),'_',num2str(36*5),'.mat'];
load(filename);
coil100 = imresize(coil100,scale);
imshow(coil100)
axes('pos',[.65 .85 .15 .15])
filename = ['coil100_',num2str(object),'_',num2str(12*5),'.mat'];
load(filename);
coil100 = imresize(coil100,scale);
imshow(coil100)
axes('pos',[.23 .02 .15 .15])
filename = ['coil100_',num2str(object),'_',num2str(48*5),'.mat'];
load(filename);
coil100 = imresize(coil100,scale);
imshow(coil100);
