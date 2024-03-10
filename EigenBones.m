function [Lambda_hat, Q, Omega, S1t, D, Dd,Thresh] = EigenBones(R,Nfft,Support,MinBoneSize)

M = size(R,1);

%--------------------------------------------------------------
%  (1) find bones
%--------------------------------------------------------------
Rf = ParaHermDFT(R,Nfft);              % DFT
D = BinwiseEVD(Rf);                 
Dd = min(-diff(D,1),[],1);             % smallest eigenvalue distance per bin;
Thresh = max(Dd)/5;                   % set threshold at 10% of max(min(EV dist.))
q = find(Dd>Thresh); 
kdiff = diff([q q(1)+Nfft]);       % '1' means index could be part of a bones (if runlength is sufficient)
BoneMarginIndices = find(kdiff>1);
BoneStart=q(1); ii=1;
for i = 1:length(BoneMarginIndices),
    BoneEnd = q(BoneMarginIndices(i)-1);
    if (BoneEnd-BoneStart)>MinBoneSize,       % check that bones aren't too short 
       Omega(ii,:) = [BoneStart BoneEnd];      % (currently this ignores that the first and last bone could be wrapped)
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
       Omega = [Omega; BoneStart BoneEnd];
   end;
end;
Q = size(Omega,1);                     % # of bones

% check if the first bone wraps around to the back
if (Omega(1,1)==1) && (Omega(Q,2)==Nfft),
   Omega(Q,2) = Omega(1,2)+Nfft;
   Omega = Omega(2:Q,:);
   Q= Q-1;
end;    

%--------------------------------------------------------------
%  (2) individual reconstruction of bones and matching
%--------------------------------------------------------------
t = (-Support:Support)';
for q = 1:Q,                               % bones should really be re-ordered --- start with longest one
   Lseg = D(:,mod( (Omega(q,1):Omega(q,2))-1,Nfft)+1);
   omega =  ( (Omega(q,1):Omega(q,2))'-1)/1024*2*pi;
   if q == 1,
      % time domain reconstruction for first bone
      S1t = inverseFourierTransformVector(Lseg',omega,t);
      S1t_norm = S1t./(ones(length(t),1)*sqrt(sum(abs(S1t).^2,1)));
   else   
      % subsequent bones
      dummy = inverseFourierTransformVector(Lseg',omega,t);
      % find permutation matrix
      dummy_norm = dummy./(ones(length(t),1)*sqrt(sum(abs(dummy).^2,1)));
      P = abs(S1t_norm'*dummy_norm);
      for i = 1:M, P(:,i)=(P(:,i)>=(max(P(:,i))-eps)); end;
      S1t = [S1t, dummy*P'];
   end;
end;

%--------------------------------------------------------------
%  (3) reconstruction --- bone-size weighted average
%--------------------------------------------------------------
W = (Omega(:,2)-Omega(:,1));
W = W/sum(W);
Lambda_hat = zeros(length(t),M);
for m = 1:M,
   Lambda_hat(:,m) = sum(S1t(:,m:M:end).*(ones(length(t),1)*W'),2);
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