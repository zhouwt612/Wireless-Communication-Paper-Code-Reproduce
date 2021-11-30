function [Pc,Pp] = RCPrecoding(H_hat,pow,t,tau2)
% Pc, Pp: precoding matrix for the common and private parts
% H_hat: the quantized channel
% tau: a factor for the imperfect CSI
% The precoding matrices are with unit power

[K,M] = size(H_hat);
Pc_orig = H_hat';
Pc = zeros(M,1);
for idx1 = 1:1:K
    Pc = Pc + Pc_orig(:,idx1);
end
Pc = sqrt(1/trace(Pc*Pc'))*Pc;
Pp_unnorm = inv(H_hat'*H_hat + (K/(pow*t*M))*((M-1+pow*t*M*tau2)/(M-1-M*tau2)*eye(M)))*H_hat';
Pp = sqrt(K/trace(Pp_unnorm*Pp_unnorm'))*Pp_unnorm;