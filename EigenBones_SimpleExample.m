% EigenBones solution to Simple Example 
% 
% Reproduces Figure 3, 4, 5 in EUSIPCO'24.
% 
% Sebastian J. Schlecht, Aalto University, 2024-03-10
% Stephan Weiss, University of Strathclyde, 2024-03-10


%--------------------------------------------------------------
%  parameters
%--------------------------------------------------------------
clear; close all;
randn('seed',0);
GeneratePEVDToyProblem3Estimates
R = Rhat2;
Nfft = 1024; M = size(R,1);
MinBoneSize = 16;
Support=3;
ThreshFraction = 1/5;

%--------------------------------------------------------------
%  EigenBones Method
%--------------------------------------------------------------

[Lambda_hat, Q, BoneBounds, BoneTime, D, Dd] = EigenBones(R,Nfft,Support,MinBoneSize,ThreshFraction);

%--------------------------------------------------------------
%  display
%--------------------------------------------------------------
%
%   plot minimum eigenvalue distance 
%
figure(1); clf;
subplot(211); plot((0:(Nfft-1))/Nfft,Dd,'b-','linewidth',1); 
hold on; 
Thresh = max(Dd)*ThreshFraction;
plot([0 1],[1 1]*Thresh,'r--','linewidth',1);
axis([0 1 0 3.2]); grid on;
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$',...
      '$5\pi/4$','$3\pi/2$','$7\pi/4$','$2\pi$'}); 
ylabel('$\delta_k$, $\cal{T}$','interpreter','latex','fontsize',10);
legend({'min.~eigenvalue distance','threshold ${\cal{T}}$'},...
       'interpreter','latex','fontsize',10,'location','NorthEast');
text(0.02, 2.5,'(a)','interpreter','latex');

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
axis([0 1 1 5.5]);
grid on;
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$',...
      '$5\pi/4$','$3\pi/2$','$7\pi/4$','$2\pi$'}); 
text(0.02, 1.5,'(b)','interpreter','latex');
xlabel('normalised angular frequency $\Omega_k$','interpreter','latex','fontsize',10);
ylabel('$\hat{\lambda}_{m,k}$',...
	'interpreter','latex','fontsize',10);
set(gcf,'OuterPosition',[230 250 570 330]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc ./figures/Segmentation.eps

%
%  plot reconstructed bones
%
figure(2);
t = (-Support:Support)';
FS = 12;
subplot(221); stem(t,real(BoneTime(:,1)),'b','linewidth',1); hold on; plot(t,imag(BoneTime(:,1)),'r*','linewidth',1);
   ylabel('$\ell_{1,q}[\tau]$','interpreter','latex','fontsize',FS);
   axis([-3.5 3.5 -1.25 3.75]); grid on;
   text(-3,2.75,'$q=1$','interpreter','latex');
subplot(223); stem(t,real(BoneTime(:,2)),'b','linewidth',1); hold on; plot(t,imag(BoneTime(:,2)),'r*','linewidth',1);
   xlabel('tag $\tau$','interpreter','latex','fontsize',FS);
   ylabel('$\ell_{2,q}[\tau]$','interpreter','latex','fontsize',FS);
   axis([-3.5 3.5 -1.25 3.75]); grid on;
   text(-3,2.75,'$q=1$','interpreter','latex');
subplot(224); 
   stem(t,real(BoneTime(:,3)),'b','linewidth',1); hold on; plot(t,imag(BoneTime(:,3)),'r*','linewidth',1);
   xlabel('tag $\tau$','interpreter','latex','fontsize',FS);
   axis([-3.5 3.5 -1.25 3.75]); grid on;
   text(-3,2.75,'$q=2$','interpreter','latex');
subplot(222); 
   % for legend only
   plot(-10,-10,'bo','linewidth',1); hold on; plot(-10,-10,'r*','linewidth',1);
   % plot actual curves
   stem(t,real(BoneTime(:,4)),'b','linewidth',1); hold on; plot(t,imag(BoneTime(:,4)),'r*','linewidth',1);
   axis([-3.5 3.5 -1.25 3.75]); grid on;
   text(-3,2.75,'$q=2$','interpreter','latex');
   legend({'$\Re\{\cdot\}$','$\Im\{\cdot\}$'},'interpreter','latex','Location','NorthEast');
set(gcf,'OuterPosition',[230 250 570 290]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc ./figures/Reconstructions.eps
   
%
%   reconstruction figure
%
figure(3); clf;
Lf = abs(fft(Lambda_hat,Nfft));
plot((0:(Nfft-1))/Nfft,Lf(:,1),'b-','linewidth',1);
hold on;
plot((0:(Nfft-1))/Nfft,Lf(:,2),'r--','linewidth',1);
axis([0 1 1 5]);
grid on;
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$',...
      '$5\pi/4$','$3\pi/2$','$7\pi/4$','$2\pi$'}); 
xlabel('normalised angular frequency $\Omega$','interpreter','latex','fontsize',10);
ylabel('$\hat{\lambda}_m(\mathrm{e}^{\mathrm{j}\Omega})$',...
	'interpreter','latex','fontsize',10);
   legend({'$m=1$','$m=2$'},'interpreter','latex','Location','NorthEast');
set(gcf,'OuterPosition',[230 250 570 220]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc ./figures/Retrieved.eps

    

