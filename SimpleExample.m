% Simple Example 
% 
% Reproduces Figure 2 in EUSIPCO'24.
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
%  eigenvalues of estimated space-time covariances
%------------------------------------------------------------------------
Lfft = 2^14;
IndexPick = [(1:Lfft/128:Lfft/8-50) (Lfft/8-49:Lfft/8+49) (Lfft/8+50:Lfft/128:Lfft*5/8-50) (Lfft*5/8-49:Lfft*5/8+49) (Lfft*5/8+50:Lfft/128:Lfft-1) Lfft];
Lhat1f = zeros(2,length(IndexPick)); Lhat2f = zeros(2,length(IndexPick));
Rhat1f = ParaHermDFT(Rhat1,Lfft); Rhat2f = ParaHermDFT(Rhat2,Lfft);
Rhat1f = Rhat1f(:,:,IndexPick); Rhat2f = Rhat2f(:,:,IndexPick);
for f = 1:length(IndexPick),
   [~,Lhat1f(:,f)] = eig(squeeze(Rhat1f(:,:,f)),'vector');
   [~,Lhat2f(:,f)] = eig(squeeze(Rhat2f(:,:,f)),'vector');
end;
Lhat1f = sort(Lhat1f,1,'descend'); Lhat2f = sort(Lhat2f,1,'descend');

%------------------------------------------------------------------------
%  display
%------------------------------------------------------------------------
FS = 12; A= 1;
Lambda1 = abs(fft(lambda1,Lfft)); Lambda2 = abs(fft(lambda2,Lfft));
%--- dummy plots for legend
h = plot([-1 -2],[-1 -2],'-','linewidth',4);
set(h(1),'color',[1 1 1]*0.75); hold on;
plot([-1 -2],[-1 -2],'b-','linewidth',A);plot([-1 -2],[-1 -2],'r-','linewidth',A);
plot([-1 -2],[-1 -2],'k--','linewidth',A);plot([-1 -2],[-1 -2],'k-','linewidth',A);
%--- now for actual curves ...
f = (0:(Lfft-1))/Lfft;
h = plot(f(IndexPick),Lambda1(IndexPick),'-','linewidth',4); 
set(h(1),'color',[1 1 1]*0.75);
hold on;
h = plot(f(IndexPick),Lambda2(IndexPick),'-','linewidth',4); 
set(h(1),'color',[1 1 1]*0.75);
plot(f(IndexPick),Lhat1f(1,:),'b--','linewidth',A); plot(f(IndexPick),Lhat1f(2,:),'r--','linewidth',A); 
plot(f(IndexPick),Lhat2f(1,:),'b-','linewidth',A); plot(f(IndexPick),Lhat2f(2,:),'r-','linewidth',A); 
axis([0 1 .5 5.5]);
xlabel('normalised angular frequency $\Omega$','interpreter','latex','fontsize',FS);
ylabel('$\lambda_{m}(\mathrm{e}^{\mathrm{j}\Omega})$, $\hat{\lambda}_{m}(\mathrm{e}^{\mathrm{j}\Omega})$',...
	'interpreter','latex','fontsize',FS);
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$','$5\pi/4$',...
      '$3\pi/2$','$7\pi/4$','$2\pi$'}); % ,...
%    'YTick',(0:1:3),'YTickLabel',...
%     {'$0$','$1$','$2$','$3$'});
grid on;
%legend({'$\lambda_m(\mathrm{e}^{\mathrm{j}\Omega})$','$\hat{\lambda}_m(\mathrm{e}^{\mathrm{j}\Omega})|_{N=10^3}$',...
%        '$\hat{\lambda}_m(\mathrm{e}^{\mathrm{j}\Omega})|_{N=10^5}$'},...
%       'interpreter','latex','fontsize',FS-2,'location','SouthEast');
legend({'ground truth','$m=1$','$m=2$','$N=10^3$','$N=10^5$'},...
       'interpreter','latex','fontsize',FS-2,'location','SouthEast');
%set(gcf,'OuterPosition',[230 250 570 280]);
%set(gca,'LooseInset',get(gca,'TightInset'));
%print -depsc Fig12a.eps
%-----  insert
Box1 = [316/512 324/512 1.5 1.65];
Box2 = [15/128 17/128 4.35 4.5];
plot(Box1([1 1 2 2 1]),Box1([3 4 4 3 3]),'k-','linewidth',A);
plot(Box2([1 1 2 2 1]),Box2([3 4 4 3 3]),'k-','linewidth',A);
axes('position',[.575 .5 .08 .4]);
   h = plot(f(IndexPick),Lambda1(IndexPick),'-','linewidth',4); 
   set(h(1),'color',[1 1 1]*0.75);
   hold on;
   h = plot(f(IndexPick),Lambda2(IndexPick),'-','linewidth',4); 
   set(h(1),'color',[1 1 1]*0.75);
   plot(f(IndexPick),Lhat1f(1,:),'b--','linewidth',A); plot(f(IndexPick),Lhat1f(2,:),'r--','linewidth',A); 
   plot(f(IndexPick),Lhat2f(1,:),'b-','linewidth',A); plot(f(IndexPick),Lhat2f(2,:),'r-','linewidth',A); 
   axis(Box1);
    set(gca,'TickLabelInterpreter','latex',...
         'XTick',[79 80 81]/128, 'XTickLabel',{'$\frac{79\pi}{128}$','$\frac{5\pi}{4}$','$\frac{81\pi}{128}$'});
     grid on;
     
axes('position',[.187 .22 .08 .35]);
   h = plot(f(IndexPick),Lambda1(IndexPick),'-','linewidth',4); 
   set(h(1),'color',[1 1 1]*0.75);
   hold on;
   h = plot(f(IndexPick),Lambda2(IndexPick),'-','linewidth',4); 
   set(h(1),'color',[1 1 1]*0.75);
   plot(f(IndexPick),Lhat1f(1,:),'b--','linewidth',A); plot(f(IndexPick),Lhat1f(2,:),'r--','linewidth',A); 
   plot(f(IndexPick),Lhat2f(1,:),'b-','linewidth',A); plot(f(IndexPick),Lhat2f(2,:),'r-','linewidth',A); 
   axis(Box2);
    set(gca,'TickLabelInterpreter','latex',...
         'XTick',[15 16 17]/128, 'XTickLabel',{'$\frac{15\pi}{128}$','$\frac{\pi}{4}$','$\frac{17\pi}{128}$'});
     grid on;
     
set(gcf,'OuterPosition',[230 250 570 350]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc ./figures/SimpleExample.eps


