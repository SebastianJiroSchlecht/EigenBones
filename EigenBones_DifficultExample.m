% EigenBones solution to Difficult Example
%
% Reproduces Figures 7, 8, 9 in EUSIPCO'24.
%
% Stephan Weiss, University of Strathclyde, 2024-03-10
% Sebastian J. Schlecht, Aalto University, 2024-03-10

%--------------------------------------------------------------
%  parameters
%--------------------------------------------------------------
clear; close all;
[R_gt,L_gt,Q_gt] = PEVDToyProblem(5);       % ground truth
% perturbation --- this must be parahermitian
randn('seed',0);
Ae = (randn(size(R_gt)) + 1i*randn(size(R_gt)))*1e-5;
R = R_gt + Ae + ParaHerm(Ae);
Nfft = 1024; M = size(R,1);
MinBoneSize = 16;
Support=3;
ThreshFraction = 1/5;

%--------------------------------------------------------------
%  EigenBones Method
%--------------------------------------------------------------

[Lambda_hat, Q, BoneBounds, BoneTime, D, Dd] = EigenBones(R,Nfft,Support,MinBoneSize, ThreshFraction);

%--------------------------------------------------------------
%  display
%--------------------------------------------------------------
%
%   plot minimum eigenvalue distance
%
figure(1); clf;
subplot(211); plot((0:(Nfft-1))/Nfft,Dd,'b-','linewidth',1);
Thresh = max(Dd)*ThreshFraction;
hold on;
plot([0 1],[1 1]*Thresh,'r--','linewidth',1);
axis([0 1 0 .6]); grid on;
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$',...
    '$5\pi/4$','$3\pi/2$','$7\pi/4$','$2\pi$'});
ylabel('$\delta_k$, $\cal{T}$','interpreter','latex','fontsize',10);
legend({'min.~eigenvalue distance','threshold ${\cal{T}}$'},...
    'interpreter','latex','fontsize',10,'location','NorthEast');

subplot(211); hold off;
[~,Lgt,~] = PEVDToyProblem(5);
Lgtf = sort(abs(PolyMatDiagSpec(Lgt,Nfft)),2,'descend');
f = (0:(Nfft-1))/Nfft;
D = sort(D,1,'descend');
h = plot(f,D(1,:),'-','linewidth',4);
set(h(1),'color',[1 1 1]*0.75);
hold on;
% only for legend
plot(-10,-10,'b-','linewidth',1);plot(-10,-10,'r--','linewidth',1);h = plot(-10,-10,'-.','linewidth',1); set(h(1),'color',[0 1 0]*0.5);
for i = 2:3,
    h = plot(f,D(i,:),'-','linewidth',4);
    set(h(1),'color',[1 1 1]*0.75);
end;
plot(f,Lgtf(:,1),'b-','linewidth',1);
plot(f,Lgtf(:,2),'r--','linewidth',1);
h = plot(f,Lgtf(:,3),'-.','linewidth',1);
set(h(1),'color',[0 1 0]*0.5);
axis([0 1 -.1 1.6]); grid on;
text(0.02, 1.3,'(a)','interpreter','latex');
ylabel('$\hat{\lambda}_{m,k}$','interpreter','latex','fontsize',10);


%
%   plot bones
%
subplot(212); hold on;
for i = 1:Q
    omega = ((BoneBounds(i,1):BoneBounds(i,2))-1);
    omega = mod(omega-1,Nfft)+1; % wrap around
    Domega = D(:,omega);
    omega(omega == Nfft) = nan; % break line at wrap around
    plot(omega/Nfft,Domega,'b-','linewidth',1);
    plot(omega([1 end])/Nfft.*ones(M,1),Domega(:,[1 end]),'b.','linewidth',5);
end
axis([0 1 -.1 1.6]);
grid on;
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$',...
    '$5\pi/4$','$3\pi/2$','$7\pi/4$','$2\pi$'});
text(0.02, 1.3,'(b)','interpreter','latex');
xlabel('normalised angular frequency $\Omega_k$','interpreter','latex','fontsize',10);
ylabel('$\hat{\lambda}_{m,k}$',...
    'interpreter','latex','fontsize',10);
set(gcf,'OuterPosition',[230 250 570 300]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc ./figures/DifficultSegmentation.eps

%
%   reconstruction figure
%
figure(3); clf;
Lf = abs(fft(Lambda_hat,Nfft));
% for legend only
plot((0:(Nfft-1))/Nfft,Lf(:,1)-100,'b-','linewidth',1);
hold on;
plot((0:(Nfft-1))/Nfft,Lf(:,2)-100,'r--','linewidth',1);
h = plot((0:(Nfft-1))/Nfft,Lf(:,3)-100,'-.','linewidth',1);
set(h(1),'color',[0 1 0]*0.5);
axis([0 1 -.1 1.6]);
h = plot(f,D(1,:),'-','linewidth',4);
set(h(1),'color',[1 1 1]*0.75);
% only for legend
plot(-10,-10,'b-','linewidth',1);plot(-10,-10,'r--','linewidth',1);h = plot(-10,-10,'-.','linewidth',1); set(h(1),'color',[0 1 0]*0.5);
for i = 2:3,
    h = plot(f,D(i,:),'-','linewidth',4);
    set(h(1),'color',[1 1 1]*0.75);
end;
Lf = abs(fft(Lambda_hat,Nfft));
plot((0:(Nfft-1))/Nfft,Lf(:,1),'b-','linewidth',1);
plot((0:(Nfft-1))/Nfft,Lf(:,2),'r--','linewidth',1);
h = plot((0:(Nfft-1))/Nfft,Lf(:,3),'-.','linewidth',1);
set(h(1),'color',[0 1 0]*0.5);
axis([0 1 -.1 1.6]);
grid on;
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$',...
    '$5\pi/4$','$3\pi/2$','$7\pi/4$','$2\pi$'});
xlabel('normalised angular frequency $\Omega$','interpreter','latex','fontsize',10);
ylabel('$\hat{\lambda}_m(\mathrm{e}^{\mathrm{j}\Omega})$',...
    'interpreter','latex','fontsize',10);
legend({'$m=1$','$m=2$','$m=3$'},'interpreter','latex','Location','NorthEast');
set(gcf,'OuterPosition',[230 250 570 230]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc ./figures/DifficultRetrieved.eps



%
%  plot reconstructed bones
%
figure(2); clf;
t = (-Support:Support)';
FS = 12; M = 3; Q = 5;
CCode = [0 0 1; 1 0 0; 0 .5 0];
seq = [1 1 3 1 1 2 3 1 3 2 3 2 2 2 3];
i = 0;
for m = 1:M,
    for q = 1:Q,
        i = i + 1;
        subplot(M,Q,(m-1)*Q+q);
        h = plot([-3.4 -3.4 3.4 3.4 -3.4],[1.15 -.45 -.45 1.15 1.15],'-','linewidth',2);
        set(h(1),'color',CCode(seq(i),:)); hold on;
        stem(t,real(BoneTime(:,m+(q-1)*M)),'b','linewidth',1); hold on;
        plot(t,imag(BoneTime(:,1)),'r*','linewidth',1);
        axis([-3.5 3.5 -.5 1.2]); grid on;
    end;
end;
subplot(M,Q,1); title('$q=1$'); ylabel('$\ell_{1,q}[\tau]$','interpreter','latex');
subplot(M,Q,2); title('$q=2$'); subplot(M,Q,3); title('$q=3$'); subplot(M,Q,4); title('$q=4$'); subplot(M,Q,5); title('$q=5$');
subplot(M,Q,6); ylabel('$\ell_{2,q}[\tau]$');
subplot(M,Q,11); ylabel('$\ell_{3,q}[\tau]$'); xlabel('$\tau$','interpreter','latex');
for i = 12:15, subplot(M,Q,i); xlabel('$\tau$','interpreter','latex'); end;
set(gcf,'OuterPosition',[230 250 570 320]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc ./figures/DifficultReconstructions2.eps

