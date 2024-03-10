% State-of-the-art Decomposition of the Simple Example
% 
% Reproduces Figure 6 in EUSIPCO'24.   
% 
% Stephan Weiss, University of Strathclyde, 2024-03-10

clear all; close all;
randn('seed',0);

%------------------------------------------------------------------------
%  ground truth parameters
%------------------------------------------------------------------------
j = sqrt(-1);
f1 = [1.6180    0.61801];  f2 = [1.6180    0.61801*j];
lambda1 = conv(f1,conj(fliplr(f1))); lambda2 = conv(f2,conj(fliplr(f2)));

%------------------------------------------------------------------------
%  estimation from data
%------------------------------------------------------------------------
Lu=1e5;
u = (randn(2,Lu)+j*randn(2,Lu))/sqrt(2);
x = zeros(2,Lu);
X(1,:) = filter(f1,1,u(1,:)); X(2,:) = filter(f2,1,u(2,:));
Rhat1 = SpaceTimeCovMatEst(X(:,1:Lu/100),1);        % 100 sample estimate
Rhat2 = SpaceTimeCovMatEst(X(:,1:Lu),1);        % 10000 sample estimate

%------------------------------------------------------------------------
%  PEVD algorithms
%------------------------------------------------------------------------
[~,Gamma1] = SMD(Rhat1,200,1e-8,1e-8);
[~,Gamma2] = SMD(Rhat2,200,1e-8,1e-8);
PSD1 = PolyMatDiagSpec(Gamma1,512);
PSD2 = PolyMatDiagSpec(Gamma2,512);
Lambda1 = PolyMatAnalyticEigValues(Rhat1);
Lambda1f = abs(fft(Lambda1.',1024));
Lambda2 = PolyMatAnalyticEigValues(Rhat2);
Lambda2f = abs(fft(Lambda2.',1024));

L1SMD = size(Gamma1,3); L2SMD = size(Gamma2,3);
L1PHEVD = size(Lambda1,2); L2PHEVD = size(Lambda2,2);
Gam1 = zeros(L1SMD,2); Gam2 = zeros(L2SMD,2);
Gam1(:,1) = squeeze(Gamma1(1,1,:)); Gam1(:,2) = squeeze(Gamma1(2,2,:));
Gam2(:,1) = squeeze(Gamma2(1,1,:)); Gam2(:,2) = squeeze(Gamma2(2,2,:));

%------------------------------------------------------------------------
%  display
%------------------------------------------------------------------------
figure(1); clf; A = 1;
subplot(221); plot((0:511)/512,abs(PSD1(:,1)),'b-','linewidth',A); 
    hold on;
    plot((0:511)/512,abs(PSD1(:,2)),'r-','linewidth',A);
    axis([0 1 1 5]); grid on;
    set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/4:1),'XTickLabel',{'$0$','$\pi/2$','$\pi$',...
      '$3\pi/2$','$2\pi$'}); 
    xlabel('norm.~ang.~freq.~$\Omega$','interpreter','latex','fontsize',10);
    ylabel('$\hat{\lambda}_m(\mathrm{e}^{\mathrm{j}\Omega})$',...
	'interpreter','latex','fontsize',10);
subplot(222); plot((0:1023)/1024,abs(Lambda1f(:,1)),'b-','linewidth',A); 
    hold on; 
    plot((0:1023)/1024,abs(Lambda1f(:,2)),'r-','linewidth',A); 
    axis([0 1 1 5]); grid on;
    set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/4:1),'XTickLabel',{'$0$','$\pi/2$','$\pi$',...
      '$3\pi/2$','$2\pi$'}); 
    xlabel('norm.~ang.~freq.~$\Omega$','interpreter','latex','fontsize',10);
    ylabel('$\hat{\lambda}_m(\mathrm{e}^{\mathrm{j}\Omega})$',...
	'interpreter','latex','fontsize',10);
subplot(223); plot(-(L1SMD-1)/2:(L1SMD-1)/2,20*log10(abs(Gam1(:,1))),'b-','linewidth',A); 
    hold on;
    plot(-(L1SMD-1)/2:(L1SMD-1)/2,20*log10(abs(Gam1(:,2))),'r--','linewidth',A);
    axis([-250 250 -120 20]); grid on;
    xlabel('lag $\tau$','interpreter','latex','fontsize',10); 
    ylabel('$\hat{\lambda}_m[\tau]$','interpreter','latex','fontsize',10);
subplot(224); 
    plot(-(L1PHEVD-1)/2:(L1PHEVD-1)/2,20*log10(abs(Lambda1(1,:))),'b-','linewidth',A); 
    hold on; 
    plot(-(L1PHEVD-1)/2:(L1PHEVD-1)/2,20*log10(abs(Lambda1(2,:))),'r--','linewidth',A);
    axis([-250 250 -120 20]); grid on;
    xlabel('lag $\tau$','interpreter','latex','fontsize',10); 
    ylabel('$\hat{\lambda}_m[\tau]$','interpreter','latex','fontsize',10);
set(gcf,'OuterPosition',[230 250 570 300]);
set(gca,'LooseInset',get(gca,'TightInset'));
%print -depsc ./figures/Benchmarks.eps

figure(2); clf;
subplot(221); 
    plot((0:511)/512,abs(PSD2(:,1)),'b-','linewidth',A); 
    hold on;
    plot((0:511)/512,abs(PSD2(:,2)),'r--','linewidth',A);
    axis([0 1 1 5]); grid on;
    set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/4:1),'XTickLabel',{'$0$','$\pi/2$','$\pi$',...
      '$3\pi/2$','$2\pi$'}); 
    text(0.02,1.5,'(a)','interpreter','latex');  
    xlabel('norm.~ang.~freq.~$\Omega$','interpreter','latex','fontsize',10);
    ylabel('$\hat{\lambda}_m(\mathrm{e}^{\mathrm{j}\Omega})$',...
	'interpreter','latex','fontsize',10);
subplot(222); 
    plot((0:1023)/1024,abs(Lambda2f(:,1)),'b-','linewidth',A); 
    hold on; 
    plot((0:1023)/1024,abs(Lambda2f(:,2)),'r--','linewidth',A); 
        axis([0 1 1 5]); grid on;
    text(0.02,1.5,'(b)','interpreter','latex');  
    set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/4:1),'XTickLabel',{'$0$','$\pi/2$','$\pi$',...
      '$3\pi/2$','$2\pi$'}); 
    xlabel('norm.~ang.~freq.~$\Omega$','interpreter','latex','fontsize',10);
    ylabel('$\hat{\lambda}_m(\mathrm{e}^{\mathrm{j}\Omega})$',...
	'interpreter','latex','fontsize',10); 
subplot(223); 
    plot(-(L2SMD-1)/2:(L2SMD-1)/2,20*log10(abs(Gam2(:,1))),'b-','linewidth',A); 
    hold on;
    plot(-(L2SMD-1)/2:(L2SMD-1)/2,20*log10(abs(Gam2(:,2))),'r--','linewidth',A);
    axis([-250 250 -120 20]); grid on;
    text(-250+250*0.04,-10,'(c)','interpreter','latex');  
    xlabel('lag $\tau$','interpreter','latex','fontsize',10); 
    ylabel('$20\log_{10}|\hat{\lambda}_m[\tau]|$','interpreter','latex','fontsize',10); 
subplot(224); 
    plot(-(L2PHEVD-1)/2:(L2PHEVD-1)/2,20*log10(abs(Lambda2(1,:))),'b-','linewidth',A); 
    hold on; plot(-(L2PHEVD-1)/2:(L2PHEVD-1)/2,20*log10(abs(Lambda2(2,:))),'r--','linewidth',A);
    axis([-250 250 -120 20]); grid on;
    text(-250+250*0.04,-10,'(d)','interpreter','latex');  
    xlabel('lag $\tau$','interpreter','latex','fontsize',10); 
    ylabel('$20\log_{10}|\hat{\lambda}_m[\tau]|$','interpreter','latex','fontsize',10); 
set(gcf,'OuterPosition',[230 250 570 300]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc ./figures/Benchmarks.eps


