function precoder = RRZF(H,Phi,sigma2,P)
% Robust Regulized Zero-forcing
[K,N] = size(H);
Phiall = zeros(N);
for idx1 = 1:1:K
    Phiall(:,:) = Phiall(:,:) + Phi(:,:,idx1);
end    
PRE = H'/(H*H'+ Phiall + sigma2/P*eye(N));
b = sqrt(P/trace(PRE*PRE'));
precoder = b*PRE;