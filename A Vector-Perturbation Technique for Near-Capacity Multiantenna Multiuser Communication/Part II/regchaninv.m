function s = regchaninv(H,u,alpha,K)
% Regularizing the inverse
s = H'/(H*H'+alpha*eye(K))*u;