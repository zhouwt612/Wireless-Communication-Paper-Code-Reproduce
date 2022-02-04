function Bbar = Bbar_f(B,f,w,QC)
% B: the effective channel matrix (NK, NK, K)
% f: precoder matrix (NK,1) f = [f1^T,...,fK^T]
% Bbar: a matrix for iteration computation
% w: weight matrix (1, K)
% QC: Quantized Channel (N, K) (same format as real channel)
% Bbar: a matrix for calculating f
[N,K] = size(QC);
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