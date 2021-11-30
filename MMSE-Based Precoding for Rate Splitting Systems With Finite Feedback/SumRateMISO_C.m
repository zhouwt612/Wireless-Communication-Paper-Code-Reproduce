function Rc = SumRateMISO_C(H,Pc,Pp,pow,t)
% Rc: sum rate of the common signal
% H: [h1,...,hK]^H channel realization
% Pc: the precoding matrix for common part
% Pp: the precoding matrix for private part
% pow: transmit power
% t: fraction of the transmitted power

[K,~] = size(H);
% Calculate the common date rate
com_interf = zeros(1,K);
gamma_ck = zeros(1,K);
for idx0 = 1:1:K
    for idx1 = 1:1:K
        com_interf(idx0) = com_interf(idx0) + abs(H(idx0,:)*Pp(:,idx1))^2;
    end
end
for idx2 = 1:1:K
    gamma_ck(idx2) = (pow*(1-t)*abs(H(idx2,:)*Pc)^2)...
        /(1 + pow*t/K*com_interf(idx2));
end
Rc = log2(1+min(gamma_ck));