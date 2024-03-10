function [Lambda_hat, Q, BoneBounds, BoneTime, D, Dd] = EigenBones(R,Nfft,Support,MinBoneSize,ThreshFraction)
%EigenBones - Compute unmajorised eigenvalues of polynomial matrix
% The method is described in EUSIPCO'24. First, the eigenvalues are
% determined by frequency bin and sorted into majorised eigenvalues.
% Secondly, bones are identified. Third, the bones are transformed via a
% partial inverse DFT and associated. Fourth, the reconstructed unmajorised
% eigenvalue is assembled.
%
% Syntax:  [Lambda_hat, Q, Omega, S1t, D, Dd, Thresh] = EigenBones(R,Nfft,Support,MinBoneSize)
%
% Inputs:
%    R - Polynomial Covariance Matrix
%    Nfft - Number of FFT points
%    Support - Time-domain support of the eigenvalues [-Support:Support]
%    MinBoneSize - Minimum Bone Size in number of frequency points
%    ThreshFraction - Eigenvalue threshold fraction of the maximum
%    difference
%
% Outputs:
%    Lambda_hat - Reconstructed unmajorised eigenvalues
%    Q - Number of Bones
%    BoneBounds - Frequency index of bone bounds
%    BoneTime - Time-domain coefficients of bones
%    D - Bin-wise eigenvalues
%    Dd - Bin-wise minimum eigenvalue difference
%
%
% Sebastian J. Schlecht, Aalto University, 2024-03-10
% Stephan Weiss, University of Strathclyde, 2024-03-10

M = size(R,1); % Number of eigenvalues

%--------------------------------------------------------------
%  (1) find bones
%--------------------------------------------------------------
Rf = ParaHermDFT(R,Nfft);              % DFT
D = BinwiseEVD(Rf);                 
Dd = min(-diff(D,1),[],1);             % smallest eigenvalue distance per bin;
Thresh = max(Dd)*ThreshFraction;       % set threshold at 10% of max(min(EV dist.))
q = find(Dd>Thresh); 
kdiff = diff([q q(1)+Nfft]);       % '1' means index could be part of a bones (if runlength is sufficient)
BoneMarginIndices = find(kdiff>1);
BoneStart=q(1); ii=1;
for i = 1:length(BoneMarginIndices),
    BoneEnd = q(BoneMarginIndices(i)-1);
    if (BoneEnd-BoneStart)>MinBoneSize,       % check that bones aren't too short 
       BoneBounds(ii,:) = [BoneStart BoneEnd];      % (currently this ignores that the first and last bone could be wrapped)
       ii = ii+1;
    end;      
    if i<length(BoneMarginIndices),
       BoneStart = q(BoneMarginIndices(i)+1);
    end;    
end;    
% check if there is a last bone 
if BoneMarginIndices(i)<length(q);
   BoneStart = q(BoneMarginIndices(i)+1); 
   BoneEnd = q(end);
   if (BoneEnd-BoneStart)>MinBoneSize, 
       BoneBounds = [BoneBounds; BoneStart BoneEnd];
   end;
end;
Q = size(BoneBounds,1);                     % # of bones

% check if the first bone wraps around to the back
if (BoneBounds(1,1)==1) && (BoneBounds(Q,2)==Nfft),
   BoneBounds(Q,2) = BoneBounds(1,2)+Nfft;
   BoneBounds = BoneBounds(2:Q,:);
   Q= Q-1;
end;    

%--------------------------------------------------------------
%  (2) individual reconstruction of bones and matching
%--------------------------------------------------------------
t = (-Support:Support)';
for q = 1:Q,                               % bones should really be re-ordered --- start with longest one
   Lseg = D(:,mod( (BoneBounds(q,1):BoneBounds(q,2))-1,Nfft)+1);
   omega =  ( (BoneBounds(q,1):BoneBounds(q,2))'-1)/1024*2*pi;
   if q == 1,
      % time domain reconstruction for first bone
      BoneTime = inverseFourierTransformVector(Lseg',omega,t);
      S1t_norm = BoneTime./(ones(length(t),1)*sqrt(sum(abs(BoneTime).^2,1)));
   else   
      % subsequent bones
      dummy = inverseFourierTransformVector(Lseg',omega,t);
      % find permutation matrix
      dummy_norm = dummy./(ones(length(t),1)*sqrt(sum(abs(dummy).^2,1)));
      P = abs(S1t_norm'*dummy_norm);
      for i = 1:M, P(:,i)=(P(:,i)>=(max(P(:,i))-eps)); end;
      BoneTime = [BoneTime, dummy*P'];
   end;
end;

%--------------------------------------------------------------
%  (3) reconstruction --- bone-size weighted average
%--------------------------------------------------------------
W = (BoneBounds(:,2)-BoneBounds(:,1));
W = W/sum(W);
Lambda_hat = zeros(length(t),M);
for m = 1:M,
   Lambda_hat(:,m) = sum(BoneTime(:,m:M:end).*(ones(length(t),1)*W'),2);
end;

end


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