function Abar = Abar_f(A,f,w,QC)
% A: the effective channel matrix (NK, NK, K)
% f: precoder matrix (NK,1) f = [f1^T,...,fK^T]
% Abar: a matrix for iteration computation
% w: weight matrix (1, K)
% QC: Quantized Channel (N, K) (same format as real channel)
% Abar: a matrix for calculating f
[N,K] = size(QC);
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