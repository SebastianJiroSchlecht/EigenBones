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
%   plot segments
%
subplot(212);
plot(((Omega(1,1):Omega(1,2))-1)/Nfft,D(:,Omega(1,1):Omega(1,2)),'b-','linewidth',1);
hold on;
plot(Omega(1,1)/Nfft*ones(1,M),D(:,Omega(1,1)),'b.','linewidth',5);
plot(Omega(1,2)/Nfft*ones(1,M),D(:,Omega(1,2)),'b.','linewidth',5);
for i = 2:Q,
    plot(((Omega(i,1):Omega(i,2))-1)/Nfft,D(:,Omega(i,1):Omega(i,2)),'b-','linewidth',1);
    plot(Omega(i,1)/Nfft*ones(1,M),D(:,Omega(i,1)),'b.','linewidth',5);
    plot(Omega(i,2)/Nfft*ones(1,M),D(:,Omega(i,2)),'b.','linewidth',5);
end;
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
        S1tu = S1t;
    else
        % subsequent sequents
        dummy = inverseFourierTransformVector(Lseg',omega,t);
        % find permutation matrix
        dummy_norm = dummy./(ones(length(t),1)*sqrt(sum(abs(dummy).^2,1)));
        P = abs(S1t_norm'*dummy_norm);
        for i = 1:M, P(:,i)=(P(:,i)>=(max(P(:,i))-eps)); end;
        S1t = [S1t, dummy*P'];
        S1tu = [S1tu dummy];
    end;
end;

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
%  plot reconstructed segments
%

figure(2); clf;
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
        stem(t,real(S1tu(:,m+(q-1)*M)),'b','linewidth',1); hold on;
        plot(t,imag(S1tu(:,1)),'r*','linewidth',1);
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
