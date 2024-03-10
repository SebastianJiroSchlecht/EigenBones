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
MinRunLength = 16;
Support=3;

%--------------------------------------------------------------
%  (1) find segments
%--------------------------------------------------------------
Rf = ParaHermDFT(R,Nfft);              % DFT
D = BinwiseEVD(Rf);                 
Dd = min(-diff(D,1),[],1);             % smallest eigenvalue distance per bin;
Thresh = max(Dd)/5;                   % set threshold at 10% of max(min(EV dist.))
k = find(Dd>Thresh); 
kdiff = diff([k k(1)+Nfft]);       % '1' means index could be part of a segment (if runlength is sufficient)
SegMarginIndices = find(kdiff>1);
SegStart=k(1); ii=1;
for i = 1:length(SegMarginIndices),
    SegEnd = k(SegMarginIndices(i)-1);
    if (SegEnd-SegStart)>MinRunLength,       % check that segments aren't too short 
       Omega(ii,:) = [SegStart SegEnd];      %   (currently this ignores that the first and last segment could be wrapped)
       ii = ii+1;
    end;      
    if i<length(SegMarginIndices),
       SegStart = k(SegMarginIndices(i)+1);
    end;    
end;    
% check if there is a last segment 
if SegMarginIndices(i)<length(k);
   SegStart = k(SegMarginIndices(i)+1); 
   SegEnd = k(end);
   if (SegEnd-SegStart)>MinRunLength, 
       Omega = [Omega; SegStart SegEnd];
   end;
end;
Q = size(Omega,1);                     % # of segments

%
%   plot minimum eigenvalue distance 
%
figure(1); clf;
subplot(211); plot((0:(Nfft-1))/Nfft,Dd,'b-','linewidth',1); 
hold on; 
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
%   plot segments
%
subplot(212);
plot(((Omega(1,1):Omega(1,2))-1)/Nfft,D(:,Omega(1,1):Omega(1,2)),'b-','linewidth',1); 
hold on;
plot(Omega(1,2)/Nfft*[1 1],D(:,Omega(1,2)),'b.','linewidth',5);
for i = 2:Q,
    plot(((Omega(i,1):Omega(i,2))-1)/Nfft,D(:,Omega(i,1):Omega(i,2)),'b-','linewidth',1); 
    if i<Q,
       plot(Omega(i,1)/Nfft*[1 1],D(:,Omega(i,1)),'b.','linewidth',5);
       plot(Omega(i,2)/Nfft*[1 1],D(:,Omega(i,2)),'b.','linewidth',5);
    else
       plot(Omega(i,1)/Nfft*[1 1],D(:,Omega(i,1)),'b.','linewidth',5);
    end;
end;
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

% check if the first segment wraps around to the back
if (Omega(1,1)==1) && (Omega(Q,2)==Nfft),
   Omega(Q,2) = Omega(1,2)+Nfft;
   Omega = Omega(2:Q,:);
   Q= Q-1;
end;    

%--------------------------------------------------------------
%  (2) individual reconstruction of segments and matching
%--------------------------------------------------------------
t = (-Support:Support)';
for k = 1:Q,                               % segments should really be re-ordered --- start with longest one
   Lseg = D(:,mod( (Omega(k,1):Omega(k,2))-1,Nfft)+1);
   omega =  ( (Omega(k,1):Omega(k,2))'-1)/1024*2*pi;
   if k == 1,
      % time domain reconstruction for first segment
      S1t = inverseFourierTransformVector(Lseg',omega,t);
      S1t_norm = S1t./(ones(length(t),1)*sqrt(sum(abs(S1t).^2,1)));
   else   
      % subsequent sequents
      dummy = inverseFourierTransformVector(Lseg',omega,t);
      % find permutation matrix
      dummy_norm = dummy./(ones(length(t),1)*sqrt(sum(abs(dummy).^2,1)));
      P = abs(S1t_norm'*dummy_norm);
      for i = 1:M, P(:,i)=(P(:,i)>=(max(P(:,i))-eps)); end;
      S1t = [S1t, dummy*P'];
   end;
end;
%
%  plot reconstructed segments
%
figure(2);
FS = 12;
subplot(221); stem(t,real(S1t(:,1)),'b','linewidth',1); hold on; plot(t,imag(S1t(:,1)),'r*','linewidth',1);
   ylabel('$\ell_{q,1}[\tau]$','interpreter','latex','fontsize',FS);
   axis([-3.5 3.5 -1.25 3.75]); grid on;
   text(-3,2.75,'$q=1$','interpreter','latex');
subplot(223); stem(t,real(S1t(:,2)),'b','linewidth',1); hold on; plot(t,imag(S1t(:,2)),'r*','linewidth',1);
   xlabel('tag $\tau$','interpreter','latex','fontsize',FS);
   ylabel('$\ell_{q,2}[\tau]$','interpreter','latex','fontsize',FS);
   axis([-3.5 3.5 -1.25 3.75]); grid on;
   text(-3,2.75,'$q=1$','interpreter','latex');
subplot(224); 
   stem(t,real(S1t(:,3)),'b','linewidth',1); hold on; plot(t,imag(S1t(:,3)),'r*','linewidth',1);
   xlabel('tag $\tau$','interpreter','latex','fontsize',FS);
   axis([-3.5 3.5 -1.25 3.75]); grid on;
   text(-3,2.75,'$q=2$','interpreter','latex');
subplot(222); 
   % for legend only
   plot(-10,-10,'bo','linewidth',1); hold on; plot(-10,-10,'r*','linewidth',1);
   % plot actual curves
   stem(t,real(S1t(:,4)),'b','linewidth',1); hold on; plot(t,imag(S1t(:,4)),'r*','linewidth',1);
   axis([-3.5 3.5 -1.25 3.75]); grid on;
   text(-3,2.75,'$q=2$','interpreter','latex');
   legend({'$\Re\{\cdot\}$','$\Im\{\cdot\}$'},'interpreter','latex','Location','NorthEast');
set(gcf,'OuterPosition',[230 250 570 290]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc ./figures/Reconstructions.eps
   
%--------------------------------------------------------------
%  (3) reconstruction --- segment-size weighted average
%--------------------------------------------------------------
W = (Omega(:,2)-Omega(:,1));
W = W/sum(W);
Lambda_hat = zeros(length(t),M);
for m = 1:M,
   Lambda_hat(:,m) = sum(S1t(:,m:M:end).*(ones(length(t),1)*W'),2);
end;

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

% save PEVDToyProblem3_NewMethodSol.mat Lambda_hat    
    
%--------------------------------------------------------------
%  some functions
%--------------------------------------------------------------
function D = BinwiseEVD(A)
M = size(A,1); Nfft = size(A,3);
D = zeros(M,Nfft);
for it = 1:Nfft,
    [~,dummy] = eig(A(:,:,it),'vector');
    % majorize
    D(:,it) = sort(real(dummy),'descend');
end
end


function x = inverseFourierTransformVector(X,w,t)

z = exp(1i*w(:));
zz = z.^(-t(:)');
x = zz \ X;

end
