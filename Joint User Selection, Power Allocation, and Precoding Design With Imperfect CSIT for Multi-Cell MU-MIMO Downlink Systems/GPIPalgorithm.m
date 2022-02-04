function precoder = GPIPalgorithm(QC,Phi,finit,w,sigma2,P,epsilon)
% generalized power iteration precoding
% QC: Quantized Channel (N, K) (same format as real channel)
% Phi: (N, N, K), the error covariance matrix for all users
% P: power constraint
itea_num = 50;
[N,K] = size(QC);
f = [finit,zeros(N*K,itea_num-1)];
for idx = 1:1:itea_num
    f(:,idx+1) = GPIPcomponent(QC,Phi,w,sigma2,P,f(:,idx));
    if norm(f(:,idx)-f(:,idx+1)) <= epsilon
        break
    end
end
PRE = f(:,idx+1);
precoder = reshape(PRE,N,K);