function precoder = GPIPalgorithm(QC,Phi,finit,w,sigma2,P,epsilon)
% generalized power iteration precoding
% QC: Quantized Channel (N, K) (same format as real channel)
% Phi: (N, N, K), the error covariance matrix for all users
% P: power constraint
itea_num = 50;
[K,N] = size(QC);
f = [finit,zeros(N*K,itea_num-1)];
for idx = 1:1:itea_num
    f(:,idx+1) = GPIPcomponent(QC,Phi,w,sigma2,P,f(:,idx));
    if norm(f(:,idx)-f(:,idx+1)) <= epsilon
        break
    end
end
PRE = f(:,idx+1);
precoder = reshape(PRE,N,K);


function f = GPIPcomponent(QC,Phi,w,sigma2,P,finit)
% f: a vector for iteration computation
[A,B] = EffecChanMtxAandB(QC,Phi,sigma2,P);
Abar = Abar_f(A,finit,w,QC);
Bbar = Bbar_f(B,finit,w,QC);
f = PrecoderGPIP(Abar,Bbar,finit);

function [A,B] = EffecChanMtxAandB(QC,Phi,sigma2,P)
% QC: Quantized Channel (N, K) (same format as real channel)
% Phi: (N, N), the error covariance matrix for all users
% sigma2: the sum of the inter-cell-interference and noise power
% P: power constraint
% A and B: the effective channel matrix (NK, NK, K)
[K,N] = size(QC);
Ainit = zeros(N*K,N*K,K);
for idx1 = 1:1:K
    for idx2 = 1:1:K
        Ainit((N*(idx2-1)+1):(N*idx2),(N*(idx2-1)+1):(N*idx2),idx1) = QC(idx1,:)'*QC(idx1,:) + Phi(:,:,idx1);
    end
end
A = Ainit + (sigma2/P)*eye(N*K,N*K);

MinusTerm = zeros(N*K,N*K,K);
for idx3 = 1:1:K
    MinusTerm((N*(idx3-1)+1):(N*idx3),(N*(idx3-1)+1):(N*idx3),idx3) = QC(idx3,:)'*QC(idx3,:);
end
B = A - MinusTerm;

function Abar = Abar_f(A,f,w,QC)
% A: the effective channel matrix (NK, NK, K)
% f: precoder matrix (NK,1) f = [f1^T,...,fK^T]
% Abar: a matrix for iteration computation
% w: weight matrix (1, K)
% QC: Quantized Channel (N, K) (same format as real channel)
% Abar: a matrix for calculating f
[K,N] = size(QC);
Aknoti = ones(1,K);
for idx1 = 1:1:K
    for idx2 = 1:1:K
        if idx2 ~= idx1
            Aknoti(1,idx1) = Aknoti(1,idx1)* (f'*A(:,:,idx2)*f);
        end
    end
end

Abarinit = zeros(N*K,N*K);
for idx3 = 1:1:K
    Abarinit(:,:) = Abarinit(:,:) + w(idx3)*(f'*A(:,:,idx3)*f)^(w(idx3)-1)*Aknoti(idx3)*A(:,:,idx3);
end
Abar = Abarinit;

function Bbar = Bbar_f(B,f,w,QC)
% B: the effective channel matrix (NK, NK, K)
% f: precoder matrix (NK,1) f = [f1^T,...,fK^T]
% Bbar: a matrix for iteration computation
% w: weight matrix (1, K)
% QC: Quantized Channel (N, K) (same format as real channel)
% Bbar: a matrix for calculating f
[K,N] = size(QC);
Bknoti = ones(1,K);
for idx1 = 1:1:K
    for idx2 = 1:1:K
        if idx2 ~= idx1
            Bknoti(1,idx1) = Bknoti(1,idx1) * (f'*B(:,:,idx2)*f);
        end
    end
end

Bbarinit = zeros(N*K,N*K);
for idx3 = 1:1:K
    Bbarinit(:,:) = Bbarinit(:,:) + w(idx3)*(f'*B(:,:,idx3)*f)^(w(idx3)-1)*Bknoti(idx3)*B(:,:,idx3);
end
Bbar = Bbarinit;

function f = PrecoderGPIP(Abar,Bbar,finit)
% f: (NK, 1) computed precoder through GPIP algorithm
% Abar: (NK, NK) a matrix for computing f
% Bbar: (NK, NK) a matrix for computing f
% finit: (NK, 1) initial f
fb = inv(Bbar)*Abar*finit;
f = fb/norm(fb,2);